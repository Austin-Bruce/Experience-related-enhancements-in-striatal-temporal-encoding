
condition = [neuronDB.condition]; 
early = find(condition==1); 
late = find(condition==2); 

intStart = -5;
intEnd = 22;
counter=1; clear PETH*; clear zPETH* 
for i_neuron = 1:length(neuronDB) % 1     
     sLE(counter) = length((neuronDB(i_neuron).events.longErrorTrialInit)); %need enough error trials, we picked > 5    
    counter=counter+1; 
end

figure(5); set(gcf, 'Color', 'white'); 
colormap('jet'); 

vector =  find(sLE>5); % more than 5 error trials
condition = [neuronDB.condition]; 

fprintf('\n Neurons with >5 errors; %d MSNs Early Training, %d MSNs Experienced',...
length(intersect(find(sLE>5), find(condition==1))), length(intersect(find(sLE>5), find(condition==2))))
fprintf('\n'); 

for i_neuron = 1:length(neuronDB)
    neuronDB(i_neuron).periEventSpikes = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.switchTrialInit, interval);
    neuronDB(i_neuron).periEventSpikesLE = peSpike(neuronDB(i_neuron).spikeTS,  neuronDB(i_neuron).events.longErrorTrialInit, interval);
end


%GLM TABLE GOES HERE 
%consider switching to 500 MS 
clear t* T* lm* pE* errorMod* ncount* a b c d e f spike* 
binSize = 0.2;
timep = (intStart+binSize:binSize:intEnd)';
tnum = sum(cellfun(@(x) size(x,1),{neuronDB.periEventSpikes})) * numel(timep); % time x trial
tnumle = sum(cellfun(@(x) size(x,1),{neuronDB.periEventSpikesLE})) * numel(timep); % time x trial
t_firing_rate = zeros(tnum,1); t_times = zeros(tnum,1); t_neuron_num = zeros(tnum,1); t_correrr = zeros(tnum,1);
t_firing_rate_le = zeros(tnumle,1); t_times_le = zeros(tnumle,1); t_neuron_num_le = zeros(tnumle,1); t_correrr_le = zeros(tnumle,1);

ncount = 1; ncountle=1; 

for i_neuron = 1:length(neuronDB)
    if sLE(i_neuron)>5;
    spike_trial = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.switchTrialInit, timep');
    spike_trial_longError = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.longErrorTrialInit, timep');
    interval_bins = intStart:binSize:intEnd;
    time_x = timep;
    
    a = histc(spike_trial',interval_bins); % time x trial
    b = a(1:end-1,:);
    nelb = ncount+numel(b)-1;

    e = histc(spike_trial_longError',interval_bins); % time x trial
    f = e(1:end-1,:);
    nelf = ncountle+numel(f)-1;

    
    t_firing_rate(ncount:nelb) = reshape(b,[],1);
    t_times(ncount:nelb) = repmat(time_x,size(b,2),1);
    t_neuron_num(ncount:nelb) = repmat(i_neuron,numel(b),1);
    t_correrr(ncount:nelb) = repmat(0,numel(b),1);
    ncount = ncount+numel(b);
    
    t_firing_rate_le(ncountle:nelf) = reshape(f,[],1);
    t_times_le(ncountle:nelf) = repmat(time_x,size(f,2),1);
    t_neuron_num_le(ncountle:nelf) = repmat(i_neuron,numel(f),1);
    t_correrr_le(ncountle:nelf) = repmat(1,numel(f),1);
    ncountle = ncountle+numel(f);
    end
end


%create a large table for the overall model
T_SPK = table(...
    [t_firing_rate; t_firing_rate_le],...
    [t_times; t_times_le],...
    [t_correrr; t_correrr_le],...
    [t_neuron_num; t_neuron_num_le],...
    'VariableNames',{'FiringRate','Times','Errors', 'Neurons'});


%create neuron by neuron glms explaining firing rate
for i_neuron = 1:length(neuronDB)
    if(sLE(i_neuron)>5)
    T_neuron = T_SPK(T_SPK.Neurons==i_neuron,:);
    lmNeuron = fitglme(T_neuron, 'FiringRate~Times+Errors');
    lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron);
    lmNeuronAnova= anova(lmNeuron); 
    pValuesRampingErrorModel(i_neuron,:) = [lmNeuronAnova{2,5}]; 
    pError(i_neuron,:) = [lmNeuronAnova{3,5}]; 
    end
end
 
clear adj* errorMod*
% Have to set FDR threshold here
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pError(find(pError>0)));  %some neurons didn't have enough errors to analyze; these are pError==0
condition = [neuronDB.condition]; 
condition=condition(pError>0);
% Chi trend for less error modulation!
errorMod = condition([find(adj_p<0.05)]');

bar([length(find(errorMod==1))./length(find(condition==1))  length(find(errorMod==2))./length(find(condition==2))])
box off; 
ylabel('Error-Modulated Fraction'); 
set(gca, 'xtick', [1:2], 'xticklabel', {'Early Training'; 'Experienced'}); 
[p,X] = chi2test([length(find(errorMod==1)) length(find(condition==1))-length(find(errorMod==1));...
    length(find(errorMod==2)) length(find(condition==2))-length(find(errorMod==2))]);
fprintf('Early %0.2f vs late %0.3f chi-square test (FDR-adjusted) = X %2.2f, p %2.2f',...
length(find(errorMod==1))./length(find(condition==1)), length(find(errorMod==2))/length(find(condition==2)), X, p); 
fprintf('\n'); 

errorCounts = [length(find(errorMod==1)) length(find(condition==1))-length(find(errorMod==1));...
    length(find(errorMod==2)) length(find(condition==2))-length(find(errorMod==2))];

csvwrite('errorCounts.csv', errorCounts');

interval = -2.5:0.25:20.5; 
i_neuron = 99; %  name: 'ER4_SPK08b_2020-04-09'
figure(51); clf; set(gcf, 'Color', 'white'); 
periEventSpike = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.switchTrialInit, interval);
periEventSpikeLongError = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.longErrorTrialInit, interval);    
periEventPlot({periEventSpike, periEventSpikeLongError}, interval, params);    % Actually plots the PE raster       

nonMod = intersect(find(pError>0.2),find([neuronDB.pTime]<0.05));


i_neuron = 48;
figure(52); clf; set(gcf, 'Color', 'white'); 
periEventSpike = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.switchTrialInit, interval);
periEventSpikeLongError = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.longErrorTrialInit, interval);    
periEventPlot({periEventSpike, periEventSpikeLongError}, interval, params);    % Actually plots the PE raster       


