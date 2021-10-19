%meta for loops 
counter = 1;
spikeDat = {};
intStart =  -5;
intEnd = 23;

Ramping = find(significance_mat_pie(:,1)==1); 
earlyRamping = intersect(Ramping, find(condiDat==1));  lengthEarlyRamping=length(earlyRamping);
lateRamping = intersect(Ramping, find(condiDat==2));

Ramping = [earlyRamping; lateRamping]; % this way, 1:lengthEarlyRamping will give you the ramping ensemble

counter=1; clear PETH; 
for i_neuron = 1:length(neuronDB) % 1
    interval = intStart:.1:intEnd; 
    periEventSpike = peSpike(neuronDB(i_neuron).spikeTS, neuronDB(i_neuron).events.switchTrialInit, interval);
    PETH(counter, :) = gksmooth(periEventSpike, interval, 0.7);  % Smooth for PCA

    counter=counter+1; 
    
    
end



%PRINCIPAL COMPONENTS ANALYSIS HERE 
cPETH = PETH(Ramping, 51:end-50); cInterval=interval(51:end-50);  
zPETH = zscore(cPETH')'; 
warning off; 
[COEFF, SCORE, ExperiencedNT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3); 
warning on; 

figure(4); clf;  set(gcf, 'Color', 'White'); 
subplot(3,1,1); 
PETHtoPlot = cPETH(1:lengthEarlyRamping,:); 
[~,sortKey] = sort(SCORE(1:lengthEarlyRamping, 1));  % For sorting by PCA
imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Early Ramping Neurons (#)'); 
xlabel('Time from Start (Seconds)'); 
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]); 
colorbar;
colormap('jet'); 

subplot(3,1,2); 
PETHtoPlot = cPETH([lengthEarlyRamping+1:end],:); 
[~,sortKey] = sort(SCORE([lengthEarlyRamping+1:end], 1)); 
imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Late Ramping Neurons (#)'); 
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]); 
xlabel('Time from Start (Seconds)'); 
colorbar;
colormap('jet'); 

subplot(3,1,3); hold on; 
plot(cInterval, COEFF(:,1)+ 0.3, 'k') ;  
plot(cInterval, COEFF(:,2)+ 0, 'k') ;  
text(16,0.4, [num2str(round(100*EXPLAINED(1)./sum(EXPLAINED))) '%'])
text(16,0.1, [num2str(round(100*EXPLAINED(2)./sum(EXPLAINED))) '%'])
xlabel('Time from Start (Seconds)'); box off; 
set(gca, 'xtick', [0 6 18], 'ytick', [0 0.3],'yticklabel', {'PC2'; 'PC1'} );

[h,p,ci,stats] = ttest2(SCORE(1:lengthEarlyRamping,1), SCORE([lengthEarlyRamping+1]:end,1));
%% 1.3
% from http://pcl.missouri.edu/sites/default/bf/two-sample.php
% N1 = 22
% N2 = 35
% t = -1.33
% Scale r = 0.707
% 
% Bayes factor in favor of the null:
% 
% Scaled JZS Bayes Factor = 1.76089
% Scaled-Information Bayes Factor = 1.288614

[h,p,ci,stats] = ttest2(SCORE(1:lengthEarlyRamping,2), SCORE([lengthEarlyRamping+1]:end,2));
%% 2.8
% from http://pcl.missouri.edu/sites/default/bf/two-sample.php
% N1 = 22
% N2 = 35
% t = -1.33
% Scale r = 0.707
% 
% Bayes factor in favor of the null:
% 
% Scaled JZS Bayes Factor = 1.76089
% Scaled-Information Bayes Factor = 1.288614



fprintf('\n |PC1| Early %0.1f(%0.1f-%0.1f)  vs Experienced: %0.1f(%0.1f-%0.1f), signrank p= %0.2f, BF 1.3  \n',...
    median(SCORE(1:lengthEarlyRamping,1)), quantile(SCORE(1:lengthEarlyRamping,1),0.25),quantile(SCORE(1:lengthEarlyRamping,1),0.75),...
    median(SCORE([lengthEarlyRamping+1]:end,1)), quantile(SCORE([lengthEarlyRamping+1]:end,1),0.25),quantile(SCORE([lengthEarlyRamping+1]:end,1),0.75),...
    ranksum(SCORE(1:lengthEarlyRamping,1), SCORE([lengthEarlyRamping+1]:end,1)));

fprintf('\n |PC2| Early %0.1f(%0.1f-%0.1f)  vs Experienced : %0.1f(%0.1f-%0.1f), signrank p= %0.2f, BF 2.8  \n',...
    median(SCORE(1:lengthEarlyRamping,2)), quantile(SCORE(1:lengthEarlyRamping,2),0.25),quantile(SCORE(1:lengthEarlyRamping,2),0.75),...
    median(SCORE([lengthEarlyRamping+1]:end,2)), quantile(SCORE([lengthEarlyRamping+1]:end,2),0.25),quantile(SCORE([lengthEarlyRamping+1]:end,2),0.75),...
    ranksum(SCORE(1:lengthEarlyRamping,2), SCORE([lengthEarlyRamping+1]:end,2)));

% figure(41)
% % To plot the PCs independently -- removed from final paper
% 
% cPETH = PETH(Ramping, 51:end-50); cInterval=interval(51:end-50);  
% zPETH = zscore(cPETH')'; 
% warning off; 
% [COEFF, SCORE, ExperiencedNT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3); 
% 
% subplot(2,2,1); 
% PETHtoPlot = cPETH(1:lengthEarlyRamping,:); 
% [~,sortKey] = sort(SCORE(1:lengthEarlyRamping, 1));  % For sorting by PCA
% imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Early Ramping Neurons (#)'); 
% xlabel('Time from Start (Seconds)'); 
% set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]); 
% colorbar;
% colormap('jet'); 
% 
% subplot(2,2,3); 
% PETHtoPlot = cPETH([lengthEarlyRamping+1:end],:); 
% [~,sortKey] = sort(SCORE([lengthEarlyRamping+1:end], 1)); 
% imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Late Ramping Neurons (#)'); 
% set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]); 
% xlabel('Time from Start (Seconds)'); 
% colorbar;
% colormap('jet'); 
%  
% 
% zPETH=cPETH(1:lengthEarlyRamping,:); 
% zPETH = zscore(zPETH')'; 
% [COEFFEarly, SCORE, ExperiencedNT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3); 
% subplot(2,2,2);
% plot(COEFFEarly(:,1:2))
% 
% % Matt wanted PCs independently
% zPETH=cPETH([lengthEarlyRamping+1:end],:);  
% zPETH = zscore(zPETH')'; 
% [COEFFLate, SCORE, ExperiencedNT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3); 
% subplot(2,2,4);
% plot(COEFFLate(:,1:2))