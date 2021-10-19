%function to identify putative medium spiny neurons and fast-spiking 
%interneurons on the basis of Berke, 2011. 
%https://www.frontiersin.org/articles/10.3389/fnsys.2011.00045/full. 
%Code written by Youngcho Kim and Austin Bruce. 

function [neuronTypes] = MSN_WF_clustering_austin(pl2_files)
 
neurons = struct;
neuronNum = 1;
for i_file = 1:length(pl2_files)
    pl2File = fullfile(pl2_files(i_file).folder,pl2_files(i_file).name);
    pl2Info = PL2GetFileIndex(pl2File);
    numSpikeChannels = numel(pl2Info.SpikeChannels);
    sortedNeuronNum =[];
    for i=1:numSpikeChannels
        tscounts = pl2Info.SpikeChannels{i}.UnitCounts;
        sortedNeuronNum(i) = sum(tscounts(2:end)>0);
    end

    % Load spike time stamps
    sortedNumCh = find(sortedNeuronNum>0);
    unitCharacters={'a','b','c','d','e','f'};
    for channelIdx = 1:length(sortedNumCh)
        for neuronIdx = 1:sortedNeuronNum(sortedNumCh(channelIdx))
            wave = internalPL2Waves(pl2File,pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name , neuronIdx);
            waveformMean = mean(wave.Waves,1);
            neuronName = sprintf('%s%s',pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name,unitCharacters{neuronIdx});
            neurons(neuronNum).fn = pl2_files(i_file).name;
            neurons(neuronNum).name = neuronName;
            neurons(neuronNum).wf.mean = waveformMean;
            neuronNum = neuronNum+1;
        end
    end
end

%% up-scaling and extract features

upsFactor = 6;
sampleRate = 60000;
timeUnitConv = 1000; % sec to milli second    
for i_neuron = 1:length(neurons)
    % waveform analysis
    waveformMean = neurons(i_neuron).wf.mean ;

    upsWf = interp(waveformMean,upsFactor); % up sampling to 40*4 KHz
    [peakAmp,peakTime] = min(upsWf);
    [hmaxAmp1, hmaxTime1] = min(abs(upsWf(1:peakTime) - peakAmp/2));
    [hmaxAmp2, hmaxTime2] = min(abs(upsWf(peakTime:end) - peakAmp/2));
    hmaxTime2 = hmaxTime2 + peakTime - 1;
    halfPeakWidthMs = abs(hmaxTime2 - hmaxTime1)/sampleRate/upsFactor*timeUnitConv;

    [troughAmp, troughTimeIdx] = max(upsWf(peakTime:end));
    troughTime = troughTimeIdx + peakTime - 1;
    peak2troughDur = (troughTime - peakTime) /sampleRate/upsFactor*timeUnitConv;
    peak2troughRatio = abs(peakAmp / troughAmp);

    wf.mean = waveformMean;
    wf.halfPeakWidthMs = halfPeakWidthMs;
    wf.peak2troughDur = peak2troughDur;
    wf.peak2troughRatio = peak2troughRatio;
    neurons(i_neuron).wf = wf;
end

%% example of waveform

i_neuron = 1;
waveformMean = neurons(i_neuron).wf.mean ;

upsWf = interp(waveformMean,upsFactor); % up sampling to 40*4 KHz
[peakAmp,peakTime] = min(upsWf);
[hmaxAmp1, hmaxTime1] = min(abs(upsWf(1:peakTime) - peakAmp/2));
[hmaxAmp2, hmaxTime2] = min(abs(upsWf(peakTime:end) - peakAmp/2));
hmaxTime2 = hmaxTime2 + peakTime - 1;
halfPeakWidthMs = abs(hmaxTime2 - hmaxTime1)/sampleRate/upsFactor*timeUnitConv;

[troughAmp, troughTimeIdx] = max(upsWf(peakTime:end));
troughTime = troughTimeIdx + peakTime - 1;
peak2troughDur = (troughTime - peakTime) /sampleRate/upsFactor*timeUnitConv;
peak2troughRatio = abs(peakAmp / troughAmp);

figure
hold on
plot(upsWf)
plot([hmaxTime2 hmaxTime1],[peakAmp/2 peakAmp/2])
text(hmaxTime1,peakAmp/2,sprintf('%s ms',num2str(halfPeakWidthMs)))
plot([peakTime troughTime],[peakAmp peakAmp])
text(peakTime,peakAmp,sprintf('%s ms',num2str(peak2troughDur)))
hold off

    
%% Clustering code

for i_neuron = 1:length(neurons)
    if ~isempty(neurons(i_neuron).wf)
        neuronWfAn(i_neuron,1) = neurons(i_neuron).wf.halfPeakWidthMs;
        neuronWfAn(i_neuron,2) = neurons(i_neuron).wf.peak2troughDur;
    end
end

% manually inspect for outlier
figure
plot(neuronWfAn(:,1),neuronWfAn(:,2),'*')


% GMM clustering
selUnits = neuronWfAn(:,1)< 0.3 & neuronWfAn(:,1) ~= 0; % exclude a two outliers

% selUnits = 1:length(neuronWfAn)% exclude a two outliers

wfFeatures = neuronWfAn(selUnits,1:2);

options = statset('Display','final','MaxIter',1500,'TolFun',1e-8);
gm = fitgmdist(wfFeatures,2,'Options',options);
idx = cluster(gm,wfFeatures);
P = posterior(gm,wfFeatures);

% 99% confidence
idxP = nan(size(idx));
idxP(P(:,1)>0.99) = 1;
idxP(P(:,2)>0.99) = 2;

% GMM uses a random seem
% peak2troughDur should be smaller for interneuron
if mean(wfFeatures(idxP==1,2)) < mean(wfFeatures(idxP==2,2)), idxInt = idxP == 2; idxNonInt = idxP == 1; idxP(idxInt) = 1; idxP(idxNonInt) = 2; end

cluster1 = (idxP == 1);
cluster2 = (idxP == 2);

figure(1);
scatter(wfFeatures(cluster1,1),wfFeatures(cluster1,2),10,P(cluster1,1),'.')
hold on
scatter(wfFeatures(cluster2,1),wfFeatures(cluster2,2),10,P(cluster2,1),'.')
hold off
legend('Cluster 1','Cluster 2')
xlim([0 0.35]); ylim([0 1]);
xlabel('halfPeakWidtH(msec)'); ylabel('peak2trough(msec)'); 
clrmap = jet(80); colormap(clrmap(9:72,:))
title('STR'); 



selIdx = find(selUnits);
for i_selIdx = 1:length(selIdx)
    if idxP(i_selIdx) == 1
        neurons(selIdx(i_selIdx)).type = 'MSN';
    elseif idxP(i_selIdx) == 2
        neurons(selIdx(i_selIdx)).type = 'INT';
    else
        neurons(selIdx(i_selIdx)).type = 'NA';
    end
end



fprintf('\nTotal MSN neurons: %d', sum(strcmp({neurons.type},'MSN')));
fprintf('\nTotal INT neurons: %d', sum(strcmp({neurons.type},'INT')));
fprintf('\nTotal Excluded neurons: %d\n', sum(strcmp({neurons.type},'NA')));

neuronTypes = neurons
% hack needed for outlier units
nonIdx = find(~selUnits); for i = 1:length(nonIdx); neuronTypes(nonIdx(i)).type = 'NA'; end;

end 
