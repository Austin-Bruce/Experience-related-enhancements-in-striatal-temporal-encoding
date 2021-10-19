%% Austin Bruce / edited by Kumar Narayanan / Checked by Eric Emmons & Youngcho Kim  
% 1/2/21
close all
clear all 

tic; 
% loadRawData_v5; %sets up dataSt and neuronDB. Raw files needed. 
LoadData %
Figure2_Revision; % 
Figure3_NeuronCounts; %Neuron by neuron generlized linear models
Figure4_PCA;
Figure5_Ensembles; %naive Bayes 
Figure6_ErrorAnalysis; %GLMe's for neurons with error-related activity
toc; 

