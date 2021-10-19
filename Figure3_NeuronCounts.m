%Utilizing GLMe's, describe proportions of neurons exhibiting time- or
%motor-dependent changes in activity.  Figure 3 in the paper. 
binSize = .200;
timep = (intStart+binSize:binSize:intEnd)';
tnum = sum(cellfun(@(x) size(x,1),{neurons.periEventSpikes})) * numel(timep); % time x trial
t_firing_rate = zeros(tnum,1); t_times = zeros(tnum,1); t_neuron_num = zeros(tnum,1); t_trialTarget = zeros(tnum, 1);
t_nosepokes_long = zeros(tnum, 1); t_nosepokes_short = zeros(tnum, 1);
ncount = 1;
trialTarget = 6;

condition=[neuronDB.condition]; 
fprintf('\n %d MSNs Naive, %d MSNs Experienced', length(find(condition==1)), length(find(condition==2))); ;

%correction for multiple comparisons calculating 
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh([[pValues.pTime], [pValues.pMotor], [pValues.pSwitch], [pValues.pShort], [pValues.pLong]]);
pVals_fdr = reshape(adj_p, [], 5); 
significance_mat = pVals_fdr < .05;
pTime_fdr = pVals_fdr(:,1); 
pMotor_fdr = pVals_fdr(:,2); 
pSwitch_fdr = pVals_fdr(:,3);
% pReward_fdr = pVals_fdr(:,4);
pShort_fdr = pVals_fdr(:,4);
pLong_fdr = pVals_fdr(:,5);


for i = 1:length(neuronDB)
    neuronDB(i).pTime = pTime_fdr(i);
    neuronDB(i).pMotor = pMotor_fdr(i);
    neuronDB(i).pSwitch = pSwitch_fdr(i);
%     neuronDB(i).pReward = pReward_fdr(i);
    neuronDB(i).pShort = pShort_fdr(i);
    neuronDB(i).pLong = pLong_fdr(i);
end 

significance_mat = pVals_fdr < .05;

switch_sig = significance_mat(:,3);
reward_sig = significance_mat(:,4);



%create pie chart describing fraction of cells with time or motor activty
%counting proportion of cells in each modulation category
significance_mat_pie = significance_mat(:,1:2);
for i = 1:length(significance_mat)
    if isequal(significance_mat_pie(i,:),[1,0])
        cellID_pie(i) = (["Ramping"]);
    end
    if isequal(significance_mat_pie(i,:),[1,1])
        cellID_pie(i) = (["Ramping + Response"]);
    end
    if isequal(significance_mat_pie(i,:),[0,1])
        cellID_pie(i) = (["Response"]);
    end
    if isequal(significance_mat_pie(i,:),[0,0])
        cellID_pie(i) = (["Nonselective"]);
    end
end

%Pie charts displaying proportion of cells selective for time or nosepoke.
ramp = [.8 .8 .8; .7 0 0;  0.8 0 0.8;0 0 0.7;];

figure(3); set(gcf, 'Color', 'white');

%find the index for the raster and plot
subplot(1,3,1); title('Ramping Neuron'); axis off; 

subplot(1,3,2); title('Response Neuron'); axis off; 




%Chi-squared test time
early_sig = significance_mat_pie(condiDat == 1)';
tm_early = sum(early_sig(:,1) == 1); %time modulated early
nm_early = sum(early_sig(:,1) == 0); %non modulated early

late_sig = significance_mat_pie(condiDat == 2)';
tm_late = sum(late_sig(:,1) == 1); %time modulated late
nm_late = sum(late_sig(:,1) == 0); %non modulated late

x1 = [repmat('a',[tm_early + nm_early],1); repmat('b',[nm_late + tm_late],1)];
x2 = [repmat(1,tm_early,1); repmat(2,nm_early,1); repmat(1,tm_late,1); repmat(2,nm_late,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2);
       
fprintf('\n Time neurons: %d (%0.2f percent) Early Training, %d (%0.2f percent)  Experienced', tm_early, tm_early/[tm_early+nm_early], tm_late, tm_late/[nm_late+tm_late]); ;
fprintf('\n Chi Square %2.1f, p %0.2f', chi2stat, pval); 

timeCounts = [tm_early nm_early; tm_late nm_late];  csvwrite('timeCounts.csv', timeCounts); 
[p,X] = chi2test(timeCounts);

%Chi-squared motor    
motor_sig = significance_mat_pie(:,2);
motor_sig_early = motor_sig(condiDat == 1);
mm_early = sum(motor_sig_early == 1); %motor modulated early
nmm_early = sum(motor_sig_early == 0); %non motor modulated early

motor_sig_late = motor_sig(condiDat == 2);
mm_late = sum(motor_sig_late == 1); %motor modulated late
nmm_late = sum(motor_sig_late == 0); %non motor modulated late

x3 = [repmat('a',[mm_early + nmm_early],1); repmat('b',[mm_late + nmm_late],1)];
x4 = [repmat(1,mm_early,1); repmat(2,nmm_early,1); repmat(1,mm_late,1); repmat(2,nmm_late,1)];
[tbl,chi2stat,pval] = crosstab(x3,x4);

fprintf('\n Motor neurons: %d (%0.2f percent) Early Training, %d (%0.2f percent)  Experienced', mm_early, mm_early/[nm_early+mm_early], mm_late, mm_late/[nmm_late+mm_late]); ;
fprintf('\n Chi Square %2.1f, p %0.2f', chi2stat, pval); 

motorCounts = [mm_early nmm_early; mm_late nmm_late];  csvwrite('motorCounts.csv', motorCounts); 
[p,X] = chi2test(motorCounts);

subplot(1,3,3); 
pH = bar(100*[tm_early/[tm_early+nm_early] tm_late/[nm_late+tm_late]; mm_early/[nm_early+mm_early] mm_late/[nm_late+mm_late]]);    
set(pH(1), 'FaceColor', [0 0.5 0]); 
set(pH(2), 'FaceColor', [0 1 0]); 
box off; 
ylabel('Modulated (%)'); 
set(gca, 'xtick', [1 2], 'xticklabel', {'Ramping'; 'Response'}); 

bothSig=sum(significance_mat_pie');
both_sig_early = bothSig(condiDat == 1);
both_early = sum(both_sig_early == 2); %time+motor modulated early
notboth_early = sum(both_sig_early < 2); %non-both  modulated early
not_either_early = sum(both_sig_early==0);

both_sig_late = bothSig(condiDat == 2);
both_late = sum(both_sig_late == 2); %time+motor modulated early
notboth_late = sum(both_sig_late < 2); %non-both modulated early
not_either_late = sum(both_sig_late==0);

x5 = [repmat('a',[both_early + notboth_early],1); repmat('b',[both_late + notboth_late],1)];
x6 = [repmat(1,both_early,1); repmat(2,notboth_early,1); repmat(1,both_late,1); repmat(2,notboth_late,1)];
[tbl,chi2stat,pval] = crosstab(x5,x6);

fprintf('\n Time+Motor neurons: %d (%0.2f percent) Early Training, %d (%0.2f percent)  Experienced', both_early, both_early/[both_early+notboth_early], both_late, both_late/[both_late+notboth_late]); ;
fprintf('\n Chi Square %2.1f, p %0.2f', chi2stat, pval)
bothCounts = [both_early notboth_early; both_late notboth_late];  csvwrite('bothCounts.csv', bothCounts); 
[p,X] = chi2test(bothCounts);


fprintf('\n Non-modulated: %d (%0.2f percent) Early Training; %d (%0.2f percent)  Experienced', not_either_early, not_either_early/sum(condiDat==1), not_either_late, not_either_late/sum(condiDat==2));
neitherCounts = [not_either_early sum(condiDat==1)-not_either_early; both_late sum(condiDat==2)-not_either_late];  csvwrite('neitherCounts.csv', bothCounts); 
[p,X] = chi2test(neitherCounts);
fprintf('\n Chi Square %2.1f, p %0.3f', X, p)


%code for generating rasters here
rasterInterval = -2:.25:20.5;

%plotting an exmample time-selective neuron
figure(31)
for i = 1:length(neuronDB)
%     if neuronDB(i).name == 'ER4_SPK07b_2020-04-09' %6-second ramp
       if neuronDB(i).name == 'EP8_SPK11b_2019-11-06' %18-second ramp
        index1 = i;
    end
end 

params.Color = [0 0 0; 1 0 0];
sT = neuronDB(index1).events.switchDeparts-neuronDB(index1).events.switchTrialInit;
[~, sortKey] = sort(sT);
periEventSpike = peSpike(neuronDB(index1).spikeTS, neuronDB(index1).events.switchTrialInit(sortKey), interval);
periEventPlot_single({periEventSpike}, rasterInterval, params);


switchNaive = sum(switch_sig(condiDat == 1));
switchExperienced = sum(switch_sig(condiDat == 2)); 
Naive = length(find(condiDat==1)); 
Experienced = length(find(condiDat==2)); 

[p,X] = chi2test([switchNaive Naive-switchNaive; switchExperienced Experienced-switchExperienced]);
fprintf('\n Switch-related neurons: %d (%0.2f percent) Early Training, %d (%0.2f percent)  Experienced, Chi2 %2.1f, p %0.2f', switchNaive, switchNaive/Naive, switchExperienced, switchExperienced/Experienced, X,p); ;
switchCounts = [switchNaive Naive-switchNaive; switchExperienced Experienced-switchExperienced];  csvwrite('switchCounts.csv', switchCounts);
csvwrite('switchCounts.csv', switchCounts');


figure(32)  
%motor-selective neuron
for i = 1:length(neuronDB)
    if neuronDB(i).name == 'EP8_SPK12a_2019-09-21'
        index1 = i;
    end
end 
params.Color = [0 0 0; 1 0 0];
periMotorSpike = peSpike(neuronDB(index1).spikeTS, sort([neuronDB(index1).events.shortPoke; neuronDB(index1).events.longPoke]), [-2:.05:2]); % Will
periEventPlot(periMotorSpike, [-2:.05:2], params);

%Wilcoxon rank-sum test to identify neurons with selectivity for short or
%long nosepokes. 
mtr_count = 1;
for i_neuron = 1:length(spikeDat)
    if pMotor_fdr(i_neuron) <.05 
        T_neuron = T_SPK(T_SPK.Neurons==i_neuron,:);
        (T_neuron.FiringRate(T_neuron.Motor == 1));
        short_fr = (T_neuron.FiringRate(T_neuron.ShortMotor == 1));
        long_fr = (T_neuron.FiringRate(T_neuron.LongMotor == 1));
        p_wilcox_shortlong{mtr_count,1} = ranksum(short_fr, long_fr);
        p_wilcox_shortlong{mtr_count,2} = neuronDB(i_neuron).name;
        p_wilcox_shortlong{mtr_count,3} = neuronDB(i_neuron).condition;
        mtr_count = mtr_count + 1;
    end
end

[h_mtr, crit_p_mtr, adj_ci_cvrg_mtr, adj_p_mtr] = fdr_bh([p_wilcox_shortlong{:,1}]);

for i = 1:length(p_wilcox_shortlong)
    p_wilcox_shortlong{i,4} = adj_p_mtr(i);
end  

side_specificity = sum([p_wilcox_shortlong{:,4}] < .05);
side_specificity_naive = sum([p_wilcox_shortlong{:,4}] < .05 & [p_wilcox_shortlong{:,3}] == 1);
side_specificity_experienced = sum([p_wilcox_shortlong{:,4}] < .05 & [p_wilcox_shortlong{:,3}] == 2);

%evaluating side-selectivity during training. 
fprintf('\n Short- or Long-Latency NP selective neurons: %d of %d motor-selective cells(%0.2f percent)', side_specificity_naive, sum([p_wilcox_shortlong{:,3}] == 1), side_specificity_naive/(sum([p_wilcox_shortlong{:,3}] == 1)));
fprintf('\n Short- or Long-Latency NP selective neurons: %d of %d motor-selective cells(%0.2f percent)', side_specificity_experienced, sum([p_wilcox_shortlong{:,3}] == 2), side_specificity_experienced/(sum([p_wilcox_shortlong{:,3}] == 2)));
sideSpecificTable = [side_specificity_naive sum([p_wilcox_shortlong{:,3}] == 1)-side_specificity_naive; side_specificity_experienced sum([p_wilcox_shortlong{:,3}] == 2)-side_specificity_experienced];
[p,X] = chi2test(sideSpecificTable); csvwrite('sideSpecificTable.csv', sideSpecificTable');
fprintf('\n Side Specific differences Chi2 %2.1f, p %0.2f', X,p); 

