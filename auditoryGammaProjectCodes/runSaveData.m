clear;close all; clc;

folderSourceString = 'N:\Projects\Divya_AuditoryProjects\data\humanECoG'; % Change this to the folder where rawData is kept

subjectName = '001PK'; expDate = '190321'; protocolName = 'GRF_001'; gridType = 'ECoG'; 

timeStartFromBaseLine = -0.6; deltaT = 2.048; 

% save data from MonkeyLogic behavior file
ML = saveMLData(subjectName,expDate,protocolName,folderSourceString,gridType);

% Read digital data from BrainProducts
[digitalTimeStamps,digitalEvents]=extractDigitalDataBrainProducts(subjectName,expDate,protocolName,folderSourceString,gridType,1.25);

% Compare ML behavior and BP files
if ~isequal(ML.allCodeNumbers,digitalEvents)
    error('Digital and ML codes do not match');
else
    clf;
    subplot(211);
    plot(1000*diff(digitalTimeStamps),'b'); hold on;
    plot(diff(ML.allCodeTimes),'r--');
    ylabel('Difference in succesive event times (ms)');
    
    subplot(212);
    plot(1000*diff(digitalTimeStamps)-diff(ML.allCodeTimes));
    ylabel('Difference in ML and BP code times (ms)');
    xlabel('Event Number');
    
    % Stimulus Onset
    stimPos = find(digitalEvents==9);
    goodStimTimes = digitalTimeStamps(stimPos);
    stimNumbers = digitalEvents(stimPos+1);
end

folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
getStimResultsML(folderExtract,stimNumbers);
goodStimNums=1:length(stimNumbers);
getDisplayCombinationsGRF(folderExtract,goodStimNums); % Generates parameterCombinations

%%%% Segement the data %%%%
getEEGDataBrainProducts(subjectName,expDate,protocolName,folderSourceString,gridType,goodStimTimes,timeStartFromBaseLine,deltaT);
