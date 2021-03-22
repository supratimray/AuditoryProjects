% This program extracts information about the auditory stimulus directly
% from ML data file
function [uniqueConditions,uniqueAudFileNames,stimData] = getAudStimInfoFromML(MLData,stimType)

if ~exist('stimType','var');            stimType=1;                     end

numTrials = length(MLData.data);

conditionNumbers = zeros(1,numTrials);
audFileNames = cell(1,numTrials);
for i=1:numTrials
    x = MLData.data(i);
    conditionNumbers(i) = x.BehavioralCodes.CodeNumbers(2);
    a = x.TaskObject.Attribute{1};
    audFileNames{i} = a{2};
end

uniqueConditions = unique(conditionNumbers);
numUniqueConditions = length(uniqueConditions);
uniqueAudFileNames = cell(1,numUniqueConditions);

for i=1:numUniqueConditions
    x = find(uniqueConditions(i)==conditionNumbers);
    
    [~,uniqueAudFileNames{i}] = fileparts(audFileNames{x(1)});
    for j=2:length(x)
        [~,tmp] = fileparts(audFileNames{x(j)});
        if ~isequal(uniqueAudFileNames{i},tmp)
            error('Condition files do not match');
        end
    end
end

if stimType==1 % AudGamma project. Find RF and RP values. All other files are listed
    typeVals = zeros(1,numUniqueConditions);
    RFVals = zeros(1,numUniqueConditions);
    RVVals = zeros(1,numUniqueConditions);
    
    for i=1:numUniqueConditions
        tmp = uniqueAudFileNames{i};
        fileDetails = strsplit(tmp,'_');
        posList = strcmpi('Type',fileDetails);
        
        if sum(posList)
            type = str2double(fileDetails{find(posList)+1});
            typeVals(i) = type;
            if type==1
                rfList = strcmpi('RF',fileDetails);
                RFVals(i) = str2double(fileDetails{find(rfList)+1});
                rvList = strcmpi('RV',fileDetails);
                RVVals(i) = str2double(fileDetails{find(rvList)+1});
            else
                disp([num2str(i) ': ' tmp ' is not a Type1 file']);
            end 
        else
            disp([num2str(i) ': ' tmp ' is not a TypeX file']);
        end
    end
    
    stimData.typeVals = typeVals;
    stimData.RFVals = RFVals;
    stimData.RVVals = RVVals;
end
end