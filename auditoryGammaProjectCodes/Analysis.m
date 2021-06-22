%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %%%%%%%%%%% ANALYSIS CODE FOR ECoG GAMMA ANALYSIS %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
[subjectName,expDate,protocolNames,good_elecs] = EcogAuditoryGammaData;
folderSourceString = 'D:\OneDrive - Indian Institute of Science\divya\NimhansRippleProject\Divya_AuditoryProjects\data\humanECoG';
gridType = 'ECoG';
timeRange = [0.25 0.75];
referenceType = 'bipolarRef'; % unipolarRef/bipolarRef/avgRef
%% Trial Sorting, Averaging, Fourier Analysis %%
TrialSortedData = cell(1,length(subjectName)); % extracted data arranged trial and electrode wise
TrialAverageData=cell(1,length(subjectName)); % data averaged across trails for each electrode
fftStim = cell(1,length(subjectName)); % Amplitude for different conditions averaged across trials
powStim = cell(1,length(subjectName)); % Power (Amplitude^2) for different conditions averaged across trials
AvgPow_acrssElecs = cell(1,length(subjectName)); % Power averaged across electrodes
powChange = cell(1,length(subjectName)); % Change in power calculated using multitaper
TimeVals  = cell(1,length(subjectName));
FreqVals  = cell(1,length(subjectName));


for ind = 1: length(subjectName)
    disp (['subjectId ' num2str(ind)])
    [TrialSortedData{ind},TrialAverageData{ind},fftStim{ind},powStim{ind},AvgPow_acrssElecs{ind},powChange{ind},TimeVals{ind},FreqVals{ind},uRf,uRv,Samp_Freq] = getValuesGamma(ind,subjectName,expDate,protocolNames,folderSourceString,gridType,referenceType,timeRange,good_elecs);
end

% first cell is subjectid, second cell  includes
% ripples(5elecs*5 frequencies*6 velocities), assr(5 elecs) and blank condition(5 elecs) (in order)
%(5th electrode is sound input)

%%
for subjectIndex = 1:length(subjectName)
    displayplots(subjectIndex,AvgPow_acrssElecs,powChange,TimeVals,FreqVals,uRf,uRv,timeRange,Samp_Freq)
end