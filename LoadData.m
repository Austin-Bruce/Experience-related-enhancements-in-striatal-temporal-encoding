%--------------------------------------------------------------------------
%called by switch_analysis_base.m 
%
%--Austin Bruce,Narayanan Lab, Sept 2020--
%--------------------------------------------------------------------------
%Description. Given .pl2 neurophysiology files and mpc behavior .txt file, 
%this script generates a) one large table of binned spike-rates around
%events of interest and b) neuron-by-neuron generalized linear models
%describing relationship between firing rate and events of interest (e.g.
%left or right motor response, "switch", reward delivery, or time-in-trial.
%--------------------------------------------------------------------------

binSize = 0.200; % in seconds

%data structures with information for files 
load('dataSt.mat')
load('neuronTypes.mat')

%"functional" vars for loops 
counter = 1;
spikeDat = {};

%choose the time window you want to look at. 0 is trial init. 18 is trial end. 
intStart = -.5;
intEnd = 18.5; %latest switch time is at 18.5 seconds, therefore the interval has been expanded +/- .5 secs. 
neuronDBcounter = 1; neuronDB=[]; 

for j_file = 1:length(dataSt)  
    file_temp.name = (dataSt(j_file).nexname(56:end));
    file_temp.name(strfind(file_temp.name, '\'))='/'; %hack for dealing with mac OS
    [file, neuronids, events] = GetFileInfo(file_temp); % This script loads the events and neuron names;
    localNeuronTypes = cell(length(neuronids),1);
    
    for i_neuronid = 1:length(neuronids)
        for i_neurontype = 1:length(neuronTypes)
            if string(neuronTypes(i_neurontype).fn) == string(dataSt(j_file).pl2_name)
                     if neuronTypes(i_neurontype).name == string(neuronids(i_neuronid)) 
                    localNeuronTypes(i_neuronid) = cellstr(neuronTypes(i_neurontype).type);
                end
            end
        end
    end
    
    %finds classes of neurons
    neuronids_wTypes = [neuronids, localNeuronTypes];
    MSNids =  neuronids_wTypes(neuronids_wTypes(:,2) == "MSN")
    FSIids =  neuronids_wTypes(neuronids_wTypes(:,2) == "INT");
    repdate = dataSt(j_file).odate;
    animal_idx = dataSt(j_file).animal_id;
    
    %extract all trial inits
    startTime = [dataSt(j_file).pl2_events.evt1.ts(1)];
    trialTypes = [dataSt(j_file).TrialAnSt.(animal_idx).programmedDuration];
    result = [dataSt(j_file).TrialAnSt.(animal_idx).outcome];
    
    % Determine if there's a "switch" response during long trials
    call_structdim = sprintf('dataSt(%d).TrialAnSt', j_file) ;
    structdim = size(getfield(eval(call_structdim), dataSt(j_file).animal_id)) ;
    loop = structdim(2);
    
    for j = 1:loop
        dataSt(j_file).TrialAnSt.(animal_idx)(j).switchYN = ~isempty(dataSt(j_file).TrialAnSt.(animal_idx)(j).SwitchDepart);
    end
    switches =  [dataSt(j_file).TrialAnSt.(animal_idx).switchYN];
 
    shortTrials = [trialTypes == 6000 & (result == 4 | result == 2)];
    longTrials =  [trialTypes == 18000 & (result == 4 | result == 2)];
    shorterrorTrials = [trialTypes == 6000 & (result == 3 | result == 1)];
    longerrorTrials = [trialTypes == 18000 & (result == 3 | result == 1)];
    switchTrials = [switches == 1];
  
    %logical loop to exlcude cases in which MedPC terminated prior to final
    %trial completion. 2 possible scenarios--animal never completed trial,
    %OR medPC protocol timed out (90minutes) during Trial. 
    if length(dataSt(j_file).TrialAnSt.(animal_idx)) == length(dataSt(j_file).pl2_events.evt3.ts);
        starts = dataSt(j_file).pl2_events.evt3.ts;
    else length(dataSt(j_file).TrialAnSt.(animal_idx)) == length(dataSt(j_file).pl2_events.evt3.ts)-1;
        starts = dataSt(j_file).pl2_events.evt3.ts(1:end-1);
    end 
    
    
    % calculate SwitchDepart and 
    animals = dataSt(j_file).animal_id;  
    trialStart = dataSt(j_file).pl2_events.evt3.ts; % CUE ON 
    trialEnd = dataSt(j_file).pl2_events.evt17.ts; % TRIAL END 
    rpInLeft = dataSt(j_file).pl2_events.evt7.ts; % LEFT RESPONSE
    rpOutLeft = dataSt(j_file).pl2_events.evt11.ts; % LEFT RELEASE
    rpInRight = dataSt(j_file).pl2_events.evt19.ts; % RIGHT RESPONSE
    rpOutRight = dataSt(j_file).pl2_events.evt15.ts;% RIGHT RELEASE
    rewards = dataSt(j_file).pl2_events.evt9.ts';% REWARD DISPENSE
    
    trialStart(find(diff(trialStart)<0.0001)+1) = [];
    
    trial = struct;
    trialNum = min(length(trialStart), length(trialEnd));
    type = typecell{strcmp(animals,typecell(:,1)),2};
    for j = 1:trialNum
        %% trial start/end time, duration, and ITI
        curTS = trialStart(j); % H 
        curTE = trialEnd(j);
        trial(j).trialStart = curTS;
        %% response time within trial
        trial(j).leftRspTimeTrial = rpInLeft(rpInLeft>curTS & rpInLeft<=curTE) - curTS;
        trial(j).leftRelTimeTrial = rpOutLeft(rpOutLeft>curTS & rpOutLeft<=curTE) - curTS;
        trial(j).rightRspTimeTrial = rpInRight(rpInRight>curTS & rpInRight<=curTE) - curTS;
        trial(j).rightRelTimeTrial = rpOutRight(rpOutRight>curTS & rpOutRight<=curTE) - curTS;

        if type == 0    % 0 indicates that a short latency trial is rewarded at the left nose poke 
            if ~isempty(trial(j).leftRelTimeTrial) && ~isempty(trial(j).rightRspTimeTrial) 
                trial(j).SwitchArrival = min(trial(j).rightRspTimeTrial(trial(j).rightRspTimeTrial > min(trial(j).leftRelTimeTrial)));
                    if ~isempty(trial(j).SwitchArrival) 
                        trial(j).SwitchDepart = max(trial(j).leftRelTimeTrial(trial(j).leftRelTimeTrial < trial(j).SwitchArrival));
                    else
                        trial(j).SwitchArrival=[];
                    end
            end
        elseif type == 1    % 1 indicates that a short latency trial is rewarded at the right nose poke             
            if ~isempty(trial(j).rightRelTimeTrial) && ~isempty(trial(j).leftRspTimeTrial) 
                trial(j).SwitchArrival = min(trial(j).leftRspTimeTrial(trial(j).leftRspTimeTrial > min(trial(j).rightRelTimeTrial)));
                    if ~isempty(trial(j).SwitchArrival) 
                        trial(j).SwitchDepart = max(trial(j).rightRelTimeTrial(trial(j).rightRelTimeTrial < trial(j).SwitchArrival));
                    else
                        trial(j).SwitchArrival=[];
                    end            
            end
        end
    end    

    shortStarts = starts(shortTrials);
    longStarts = starts(longTrials);
    switchStarts = starts(switchTrials);
    shortErrorStarts = starts(shorterrorTrials);
    longErrorStarts = starts(longerrorTrials);
    sw_trial = ~cellfun(@(x) isempty(x), {trial.SwitchDepart});
    switchDeparts = [trial(sw_trial).trialStart]' + [trial(sw_trial).SwitchDepart]';
    switchArrivals = [trial(sw_trial).trialStart]' + [trial(sw_trial).SwitchArrival]';
    neuronIdx = 1;
    neuron_count{j_file} = length(MSNids);
    
    %% Loads neural data from pl2 file; code by YC Kim / Austin
    %% Can be ~10 ms error.
    %         % Loop to load the events
    neurons = struct;
    for i_neuron = 1:length(MSNids) % can be changed if interested in FSIids
        neurons(neuronIdx).name = MSNids{i_neuron};
        [~, spike_ts] = nex_ts(file.name, MSNids{i_neuron}); % Loads the timestamps
        neurons(neuronIdx).spikeTS = spike_ts;  % Copies to neuron struct
        neurons(neuronIdx).numTS = length(spike_ts);  % Gets the number of timestamps
        neurons(neuronIdx).events.short =  shortStarts;
        neurons(neuronIdx).events.long =  longStarts;
        neurons(neuronIdx).events.switches = switchStarts;
        neurons(neuronIdx).events.switchTimes = switchDeparts;
        neurons(neuronIdx).events.switchArrives = switchArrivals;
        neurons(neuronIdx).events.rewards = rewards;

        if ~isempty(neurons(neuronIdx).events)
            interval = intStart:binSize:intEnd;
            periEventSpike = peSpike(neurons(neuronIdx).spikeTS, neurons(neuronIdx).events.switches, interval);
            periEventMotor =  peSpike([dataSt(j_file).pl2_events.evt19.ts; dataSt(j_file).pl2_events.evt7.ts], neurons(neuronIdx).events.switches, interval);
            % Need to analyze short vs. long nosepokes separately?
            periEventRewards = peSpike(dataSt(j_file).pl2_events.evt9.ts, neurons(neuronIdx).events.switches, interval);
            periEventSwitchTimes = peSpike(neurons(neuronIdx).events.switchTimes, neurons(neuronIdx).events.switches, interval);
            
            periEventRight = peSpike(dataSt(j_file).pl2_events.evt19.ts, neurons(neuronIdx).events.switches, interval);
            periEventLeft = peSpike(dataSt(j_file).pl2_events.evt7.ts, neurons(neuronIdx).events.switches, interval);
            
            
            if ~isempty(periEventSpike)
                PETH(counter, :) = gksmooth(periEventSpike, interval, .3);
            else
                PETH(counter, :) = zeros(1, length(PETH));
            end
        end
        
        neurons(neuronIdx).periEventSpikes = periEventSpike;  % Makes peri-event spike array
        
        %add relevant behavioral events
        spikeDat{counter}= periEventSpike;
        neuronNames{counter}= strcat(animal_idx,'_',neurons(neuronIdx).name, '_', dataSt(j_file).odate);
        rewardDat{counter} = periEventRewards;
        condiDat(counter) = dataSt(j_file).condition; %1 is early, 2 is late
        switchDat{counter} = periEventSwitchTimes;
        motorDat{counter} = periEventMotor;
        leftDat{counter} = periEventLeft;
        rightDat{counter} = periEventRight;
        
        neuronDB(counter).name = neuronNames{counter};
        neuronDB(counter).nexfile = file_temp.name;
        neuronDB(counter).spikeTS = spike_ts;  % Copies to neuron struct
        neuronDB(counter).numTS = length(spike_ts);  % Gets the number of timestamps
        neuronDB(counter).condition = dataSt(j_file).condition;
        neuronDB(counter).events.shortTrialInit =  shortStarts;
        neuronDB(counter).events.longTrialInit =  longStarts;
        neuronDB(counter).events.switchTrialInit = switchStarts;
        neuronDB(counter).events.switchDeparts = switchDeparts;
        neuronDB(counter).events.switchArrives = switchArrivals;
        neuronDB(counter).events.Rewards = rewards;
        neuronDB(counter).events.shortErrorTrialInit = shortErrorStarts;
        neuronDB(counter).events.longErrorTrialInit = longErrorStarts;
        
        if ismember(animal_idx, ['ER3', 'ER4', 'ER5', 'EP8'])
            neuronDB(counter).events.shortPoke = dataSt(j_file).pl2_events.evt19.ts;
            neuronDB(counter).events.longPoke = dataSt(j_file).pl2_events.evt7.ts;
            neuronDB(counter).events.shortRelease = dataSt(j_file).pl2_events.evt15.ts;
            neuronDB(counter).events.longRelease = dataSt(j_file).pl2_events.evt11.ts;            
            shortDat{counter} = rightDat{counter};
            longDat{counter} = leftDat{counter};
            
        else
            neuronDB(counter).events.shortPoke = dataSt(j_file).pl2_events.evt7.ts;
            neuronDB(counter).events.longPoke = dataSt(j_file).pl2_events.evt19.ts;
            neuronDB(counter).events.shortRelease = dataSt(j_file).pl2_events.evt11.ts;
            neuronDB(counter).events.longRelease = dataSt(j_file).pl2_events.evt15.ts;
            shortDat{counter} = leftDat{counter};
            longDat{counter} = rightDat{counter};
        end
        
        counter=counter+1;  % over whole ensemble
        neuronIdx = neuronIdx+1;  % over each animal
    end
    
end


timep = (intStart+binSize:binSize:intEnd)';
tnum = sum(cellfun(@(x) size(x,1),{neurons.periEventSpikes})) * numel(timep); % time x trial
t_firing_rate = zeros(tnum,1); t_times = zeros(tnum,1); 
t_neuron_num = zeros(tnum,1); t_trialTarget = zeros(tnum, 1);
t_nosepokes_long = zeros(tnum, 1); t_nosepokes_short = zeros(tnum, 1);
ncount = 1;
trialTarget = 6;


for i_neuron = 1:length(spikeDat)
%     spike_trial = neurons(i_neuron).periEventSpikes;
    spike_trial = double(spikeDat{i_neuron});
    short_rsps = double(shortDat{i_neuron});
    long_rsps = double(longDat{i_neuron});
    left_rsps = double(leftDat{i_neuron});
    right_rsps = double(rightDat{i_neuron});
    motor_trial = double(motorDat{i_neuron});
%     reward_trial = diag(rewardDat{i_neuron}(:,1));
    interval_bins = intStart:binSize:intEnd;
    switch_rsps = switchDat{i_neuron};
    time_x = timep;
    
    a = histc(spike_trial',interval_bins); % time x trial
    b = a(1:end-1,:);
    nelb = ncount+numel(b)-1;
    
    y = histc(motor_trial',interval_bins); %np x trial
    dimy = size(y);
    if dimy(1) > 1
        z = y(1:end-1,:);
    else 
        z = y(:,1:end-1);
        z = z';
    end
    
    %short rsps
    f = histc(short_rsps', interval_bins);
    dimf = size(f);
    if dimf(1) > 1
        g = f(1:end-1,:)';
    else
        g = f(:,1:end-1)';
        g = g';
    end
    
    %long rsps
    h = histc(long_rsps', interval_bins);
    dimh = size(h);
    if dimh(1) > 1
        i = h(1:end-1,:);
    else
        i = h(:,1:end-1);
        i = i';
    end
    
    %left rsps
    j = histc(left_rsps', interval_bins);
    dimj = size(j);
    if dimj(1) > 1
        k = j(1:end-1,:);
    else
        k = j(:,1:end-1);
        k = k';
    end
    
    %right rsps
    l = histc(right_rsps', interval_bins);
    diml = size(l);
    if diml(1) > 1
        m = l(1:end-1,:);
    else
        m = l(:,1:end-1);
        m = m';
    end
    
    
%     
%     reward_trial(reward_trial == 0) = NaN;
%     r = histc(reward_trial', interval_bins);
%     s = r(1:end-1,:);

    switch_rsps  = diag(switch_rsps);
    switch_rsps(switch_rsps== 0) = NaN;
    d = histc(switch_rsps, interval_bins);
    e = d(1:end-1,:);
    
    
    t_firing_rate(ncount:nelb) = reshape(b,[],1);
    t_times(ncount:nelb) = repmat(time_x,size(b,2),1);
    t_targets = t_times == 6;
    t_neuron_num(ncount:nelb) = repmat(i_neuron,numel(b),1);
    t_nosepokes(ncount:nelb) = reshape(z,[],1);
%     t_rewards(ncount:nelb) = reshape(s,[],1);
    t_switches(ncount:nelb) = reshape(e,[],1);
    t_short(ncount:nelb) = reshape(g,[],1);
    t_long(ncount:nelb) = reshape(i,[],1);
    t_left(ncount:nelb) = reshape(k,[],1);
    t_right(ncount:nelb) = reshape(m,[],1);
    ncount = ncount+numel(b);
end

% 

t_switches(isnan(t_switches)) = 0;
% t_rewards(isnan(t_rewards)) = 0;
if size(t_nosepokes,1)==1; t_nosepokes = t_nosepokes'; end;
% if size(t_rewards,1)==1; t_rewards = t_rewards'; end;
if size(t_switches,1)==1; t_switches = t_switches'; end;
if size(t_left,1)==1; t_left = t_left'; end;
if size(t_right,1)==1; t_right = t_right'; end;
if size(t_short,1)==1; t_short = t_short'; end;
if size(t_long,1)==1; t_long = t_long'; end;


%create a large table for the overall model
% T_SPK = table(t_firing_rate,t_times,t_neuron_num, t_nosepokes, t_switches, t_rewards, 'VariableNames',{'FiringRate','Times','Neurons', 'Motor', 'Switch', 'Reward'});
T_SPK = table(t_firing_rate,t_times,t_neuron_num, t_nosepokes, t_switches, t_short, t_long, t_left, t_right, 'VariableNames',{'FiringRate','Times','Neurons', 'Motor', 'Switch', 'ShortMotor', 'LongMotor', 'LeftMotor', 'RightMotor'});
lm = fitglme(T_SPK, 'FiringRate~Times+Motor+Switch+(1|Neurons)');
anova(lm);  

pValues = struct();
pValues.names = neuronNames;

%create neuron by neuron glme's and store the result in one struct
for i_neuron = 1:length(spikeDat)
    T_neuron = T_SPK(T_SPK.Neurons==i_neuron,:);
    
    lmNeuron = fitglme(T_neuron, 'FiringRate~Times');
    lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron); lmNeuronAnova = anova(lmNeuron); 
    pValues(i_neuron).pTime = [lmNeuronAnova{2,5}];
    
    lmNeuron = fitglme(T_neuron, 'FiringRate~Motor');
    lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron); lmNeuronAnova = anova(lmNeuron); 
    pValues(i_neuron).pMotor = [lmNeuronAnova{2,5}];
    
    lmNeuron = fitglme(T_neuron, 'FiringRate~ShortMotor');
    lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron); lmNeuronAnova = anova(lmNeuron); 
    pValues(i_neuron).pShort = [lmNeuronAnova{2,5}];

    lmNeuron = fitglme(T_neuron, 'FiringRate~LongMotor');
    lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron); lmNeuronAnova = anova(lmNeuron); 
    pValues(i_neuron).pLong = [lmNeuronAnova{2,5}];

    lmNeuron = fitglme(T_neuron, 'FiringRate~Switch');
    lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron); lmNeuronAnova = anova(lmNeuron); 
    pValues(i_neuron).pSwitch = [lmNeuronAnova{2,5}];
    
%     lmNeuron = fitglme(T_neuron, 'FiringRate~Reward');
%     lmNeuronbyNeuronAnova{i_neuron} = anova(lmNeuron); lmNeuronAnova = anova(lmNeuron); 
%     pValues(i_neuron).pReward = [lmNeuronAnova{2,5}];
end



%params for rasters
nexFiles = struct;
params = struct;
params.Raster = 'on';
params.CI = 'off';
params.Bar = 'off';
params.ZeroLine = 'on';
params.RandPerm = 0;
params.Color = [ 0 0 0; 1 0 0; 0 0 1; 1 0 1];  % don't need a colormap command