
condition = [neuronDB.condition]; 
early = find(condition==1); 
late = find(condition==2); 
earlyNeurons = neuronDB(early); 
lateNeurons=neuronDB(late);

earlyStats=EnsembleAnalysisSwitch(earlyNeurons, 0.5, 1);  % Currently, 20 trial minimum; does better with more trials, obviously
lateStats=EnsembleAnalysisSwitch(lateNeurons, 0.5, 1);  % Currently, 20 trial minimum; does better with more trials, obviously


% Reviewer wondered about main effect of bin size; 
bins = [0.1 0.2 0.5 1 2 3]; t = [];
for i_bin = 1:6
    binStats(i_bin)=EnsembleAnalysisSwitch(lateNeurons, 0.5, bins(i_bin));  % Currently, 20 trial minimum; does better with more trials, obviously
    b = repmat(bins(i_bin), 1, length(binStats(i_bin).r2));
    r =  binStats(i_bin).r2; 
    t = [b' r'; t]; 
end
t = table(t(:,1), t(:,2), 'VariableNames', {'binSize', 'R2'});
lmBins = fitglm(t, 'binSize~R2');  % no main effect
signrank( binStats(2).r2,  binStats(4).r2) % 200 ms is the same as 1 seconds




figure(5); set(gcf, 'color', 'white'); % Lots of cosmetic stuff possibly, but a quick check
subplot(3,2,1); imagesc(earlyStats.times, earlyStats.times, earlyStats.y, [0 1.0]); 
set(gca, 'xtick', [0 6 18], 'ytick', [0 6 18]); 
colorbar; ylabel('Predicted Time (Seconds)'), xlabel('Observed Time (Seconds)'); title('Naive'); 
subplot(3,2,2); imagesc(lateStats.times, lateStats.times, lateStats.y, [0 1.0]); 
set(gca, 'xtick', [0 6 18], 'ytick', [0 6 18]); 
colorbar; ylabel('Predicted Time (Seconds)'), xlabel('Observed Time (Seconds)'); title('Experienced'); 


subplot(3,2,3:4); cla; hold on; 

clear vdata;
vdata{1} = [earlyStats.r2]'; 
vdata{2} = [earlyStats.shuffr2]'; 
vdata{3} = [lateStats.r2]'; 
vdata{4} = [earlyStats.shuffr2]'; 

labels = [repmat('Early',length(vdata{1}), 1); repmat('Eshuf',length(vdata{2}), 1);  repmat('Late ',length(vdata{3}), 1); repmat('Lshuf',length(vdata{4}), 1);  ];
categories = [ones(1,length(vdata{1})) 2*ones(1,length(vdata{2})) 3*ones(1,length(vdata{2})) 4*ones(1,length(vdata{4}))];

vdata = ([vdata{1}; vdata{2}; vdata{3}; vdata{4} ]);
cmap = [ .2 1 .2; 0.1 0.3 0.1];
vfig = violinplot(vdata, categories,'ShowMean', true, 'ViolinAlpha', 0.1);
vfig(1).ViolinColor=[0.1 0.3 0.1]; 
vfig(2).ViolinColor=[0.7 0.7 0.7]; 
vfig(3).ViolinColor=[.2 1 .1]; 
vfig(4).ViolinColor=[0.7 0.7 0.7]; 
set(gca, 'xtick', 1:4, 'xticklabel', {'Early'; 'Shuff'; 'Late'; 'Shuff'}); 
ylabel('R^2'); 

fprintf('\n Stats:_________________'); 
fprintf('\n R2: Early %2.2f(%0.2f-%0.2f) Vs Late %2.2f(%0.2f-%0.2f)', median(earlyStats.r2), quantile(earlyStats.r2,0.25),quantile(earlyStats.r2,0.75), median(lateStats.r2), quantile(lateStats.r2, 0.25), quantile(lateStats.r2, 0.75))
fprintf('\n Shuff R2: Early %2.2f(%0.2f-%0.2f) Vs Late %2.2f(%0.2f-%0.2f)', median(earlyStats.shuffr2), quantile(earlyStats.shuffr2, 0.25), quantile(earlyStats.shuffr2, 0.75), median(lateStats.shuffr2), quantile(lateStats.shuffr2,0.25), quantile(lateStats.shuffr2,0.75))

fprintf('\n R2: signrank Early V. Late p = %0.3f, cohen d %2.1f', signrank(earlyStats.r2, lateStats.r2), cohend(earlyStats.r2, lateStats.r2))
fprintf('\n R2: signrank Early Vs Shuff p = %0.4f', signrank(earlyStats.r2, earlyStats.shuffr2))
fprintf('\n R2: signrank Late Vs Shuff p = %0.4f', signrank(lateStats.r2, lateStats.shuffr2))
fprintf('\n'); 
csvwrite('earlyStats.csv', earlyStats.r2')
csvwrite('lateStats.csv', lateStats.r2')
csvwrite('earlyStatsShuff.csv', earlyStats.shuffr2')
csvwrite('lateStatsShuff.csv', lateStats.shuffr2')


TimeNeurons =  neuronDB(find(cellID_pie=='Ramping'&condition==2));
MotorNeurons =  neuronDB(find(cellID_pie=='Response'&condition==2));
TimeMotorNeurons =  neuronDB(find(cellID_pie=='Ramping + Response'&condition==2));

NoTime = find(cellID_pie~='Ramping');
NoTimeMotor = find(cellID_pie~='Ramping + Response');
NoTime = intersect(NoTime, NoTimeMotor); 
% Experienced = find(condition==2); 
Experienced = find(condition==2);
NoTimeNeurons =  neuronDB(intersect(NoTime, Experienced));


rampStats=EnsembleAnalysisSwitch(TimeNeurons, 0.5, 1);  % Currently, 8 trial minimum; does better with more trials, obviously
motorNeuronsStats=EnsembleAnalysisSwitch(MotorNeurons, 0.5, 1);  % Currently, 8 trial minimum; does better with more trials, obviously
bothNeuronsStats=EnsembleAnalysisSwitch(TimeMotorNeurons, 0.5, 1);  % Currently, 8 trial minimum; does better with more trials, obviously
noTimeStats=EnsembleAnalysisSwitch(NoTimeNeurons, 0.5, 1);  % Check that I did this right...


csvwrite('rampStats.csv', rampStats.r2')
csvwrite('bothNeuronsStats.csv', bothNeuronsStats.r2')
csvwrite('motorNeuronsStats.csv', motorNeuronsStats.r2')
csvwrite('noTimeStats.csv', noTimeStats.r2')


subplot(3,2,5:6);  cla; hold on; 
clear vdata;
vdata{1} = [bothNeuronsStats.r2]'; 
vdata{2} = [rampStats.r2]'; 
vdata{3} = [motorNeuronsStats.r2]'; 
vdata{4} = [noTimeStats.r2]'; 

categories = [ones(1,length(vdata{1})) 2*ones(1,length(vdata{2})) 3*ones(1,length(vdata{3})) 4*ones(1,length(vdata{4})) ];
vdata = ([vdata{1}; vdata{2}; vdata{3}; vdata{4};  ]);
vfig = violinplot(vdata, categories,'ShowMean', true, 'ViolinAlpha', 0.1);
vfig(1).ViolinColor=[0 0.7 0]; 
vfig(2).ViolinColor=[0 0.7 0.7]; 
vfig(3).ViolinColor=[0.7 0.7 0]; 
vfig(4).ViolinColor=[0.1 0.1 0.1]; 

set(gca, 'xtick', 1:4, 'xticklabel', {'Time+Motor'; 'Time'; 'Motor'; 'No Time'}); 
ylabel('R^2'); 
set(gcf, 'Color', 'white');
fprintf('\n Stats:_________________'); 
fprintf('\n R2: BothStats %2.2f(%0.2f-%0.2f); Ramp %2.2f(%0.2f-%0.2f);  Motor %2.2f(%0.2f-%0.2f): No Ramp %2.2f(%0.2f-%0.2f)',...
median(bothNeuronsStats.r2), quantile(bothNeuronsStats.r2,0.25),quantile(bothNeuronsStats.r2,0.75),...
median(rampStats.r2), quantile(rampStats.r2,0.25),quantile(rampStats.r2,0.75),...
median(motorNeuronsStats.r2), quantile(motorNeuronsStats.r2,0.25),quantile(motorNeuronsStats.r2,0.75),...
median(noTimeStats.r2), quantile(noTimeStats.r2,0.25),quantile(noTimeStats.r2,0.75));

fprintf('\n R2: signrank Time+Motor vs Time = %0.2f,  Bayes Factor d %2.1f', signrank(rampStats.r2, bothNeuronsStats.r2), 1/BayesFactor(rampStats.r2, bothNeuronsStats.r2))
fprintf('\n'); 
fprintf('\n R2: signrank Time+Motor Vs Motor p = %0.5f,  cohen d %2.1f', signrank(bothNeuronsStats.r2, motorNeuronsStats.r2), cohend(bothNeuronsStats.r2, motorNeuronsStats.r2))
fprintf('\n R2: signrank Time Vs Motor p = %0.3f,  cohen d %2.1f', signrank(rampStats.r2, motorNeuronsStats.r2), cohend(rampStats.r2, motorNeuronsStats.r2))
fprintf('\n R2: signrank Time+Motor Vs No Time p = %0.3f,  cohen d %2.1f', signrank(bothNeuronsStats.r2, noTimeStats.r2), cohend(bothNeuronsStats.r2, noTimeStats.r2))
fprintf('\n R2: signrank Time Vs No Time p = %0.4f,  cohen d %2.1f', signrank(rampStats.r2, noTimeStats.r2), cohend(rampStats.r2, noTimeStats.r2))
fprintf('\n'); 

