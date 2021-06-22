%%%% This function arranges extracted data for further analysis %%%%%

%%%% Input arguments include 
% protocol related details - (1) index - subject for which you want to run the
% function. (2)subjectName.(3)expDate.(4)protocolName
%folderSourceString - string that gives location/path of stored extracted data
% gridType - Recording scale(EcoG/EEG)
%refernceType = Can be either one of the three - 'unipolarRef','bipolarRef','avgRef'
% timeRange - time values for which analysis is to be done - kept
% 0.25-0.75s (default - both stimulus and baseline periods)
%good_elecs - Electrode numbers for which power is to be averaged across
%trials. 
%useERP - 0(taking fft across each trial and then averaging across trials) or 1(averaging across each trial and then doing fft)

function [TrialSortedData,TrialAverageData,fftStim,powStim,avgPow_acrssElecs,powChange,TimeVals,FreqVals,uRf,uRv,Fs] = getValGamma(index,subjectName,expDate,protocolNames,folderSourceString,gridType,referenceType,timeRange,good_elecs,useERP,params)

if ~exist('useERP','var'); useERP = 0; end
if ~exist('timeRange','var');timeRange = [0.25 0.75];end

%% load essential data
folderName = string(fullfile(folderSourceString,'data',subjectName{index},gridType,expDate{index},protocolNames{index}));
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
x = load(fullfile(folderLFP,'lfpInfo.mat'));
ElecSet=sort(x.analogChannelsStored);
timeVals = x.timeVals;
p = load(fullfile(folderExtract,'parameterCombinations.mat'));
MLData = load(fullfile(folderExtract,'ML.mat'));
if exist('badTrials','var')
    load(fullfile(folderSegment,'badTrials'));%#ok<LOAD>
else
    disp ('bad trial file doesnot exists')
    badTrials=[];
end

%%
%%retrieving stimulus information from parameterCombinations%%
[uniqueConditions,~,stimData] = getAudStimInfoFromML(MLData);

if isequal(uniqueConditions,p.oValsUnique)
    numUniqueConditions = length(uniqueConditions);
    parameterCombinations = cell(1,numUniqueConditions);
    
    for ipc=1:numUniqueConditions
        %sorting trial numbers for each condition
        parameterCombinations{ipc} = p.parameterCombinations{1,1,1,1,ipc};
    end
end

uRf = unique(stimData.RFVals); %identifying unique ripple frequencies
uRv = unique(stimData.RVVals); %identifying unique ripple velocities

RfVals = stimData.RFVals(1:length(uRf)*length(uRv));%   ripple conditions
RvVals = stimData.RVVals(1:length(uRf)*length(uRv));%  ripple conditions
types = unique(stimData.typeVals) ; %type of conditions


avgAnalogData = cell(1,length(types));
fftStim= cell(1,length(types));
powStim = cell(1,length(types));
avgPow_acrssElecs = cell(1,length(types)); % power averaged across good_elecs
powChange = cell(1,length(types)-1); % change in power averaged across good_elecs (stim-base conditions)
%% running 'for' loop for each electrode %%
if  strcmpi(referenceType,'bipolarRef')
    disp('bipolarRef')
    load(fullfile(folderSegment,'LFP','BipolarRef.mat'));
    ElecSet = 1:length(bipolarData);
else
    ElecSet=sort(x.analogChannelsStored);
    if strcmpi(referenceType,'avgRef')
        disp('Average Referencing')
    else
        disp('Unipolar Referencing')
    end
end

%% sorting trial wise data and averaging it for different conditions
analogDataAllGra = cell(1,length(types)); % one cell for ripple, assr and blank each
analogDataAllGra{1}= cell(length(ElecSet),length(uRf),length(uRv)); %hardcoded - creating cell for arranging ripple data

for ielec = 1:length(ElecSet)
    disp(['elec' num2str(ElecSet(ielec))]);
    
    %%Loading data for one electrode at a time
    clear analogData
    
    %% referencing
    if strcmpi(referenceType,'unipolarRef')
        load(fullfile(folderSegment,'LFP',['elec' num2str(ElecSet(ielec))]),'analogData');
        
    elseif strcmpi(referenceType,'avgRef')
        load(fullfile(folderSegment,'LFP',['elec' num2str(ElecSet(ielec))]),'analogData');
        x = load(fullfile(folderSegment,'LFP','AvgRef.mat'));
        analogData = analogData - x.analogData;
        
    else
        
        analogData = bipolarData{ielec};
    end
    

    for type = 1%ripple conditions
        for irf = 1:length(unique(RfVals))
            for irv = 1:length(unique(RvVals))
                ind = find(RvVals == RvVals(irv) & RfVals == uRf(irf));
                %disp(['condNum ' num2str(ind)])
                trialNums = parameterCombinations{ind};
                trialNums = setdiff(trialNums,badTrials);
                analogDataAllGra{1}{ielec,irf,irv} = cat(1,analogDataAllGra{1}{ielec,irf,irv},analogData(trialNums,:));
                avgAnalogData{1}{ielec,irf,irv} = mean(analogDataAllGra{1}{ielec,irf,irv}); %averaging across trials
            end
        end
    end
    
    for type = 2 % assr
        trialNums = parameterCombinations{ind+1};
        trialNums = setdiff(trialNums,badTrials);
        analogDataAllGra{2}{ielec,1}= (analogData(trialNums,:));
        avgAnalogData{2}{ielec,1} = mean(analogDataAllGra{2}{ielec,1} );
    end
    
    for type = 0 % blank/baseline
        trialNums = parameterCombinations{ind+2};
        trialNums = setdiff(trialNums,badTrials);
        analogDataAllGra{3}{ielec,1}= (analogData(trialNums,:));
        avgAnalogData{3}{ielec,1} = mean(analogDataAllGra{3}{ielec,1});
    end
    
    %creating frequency axis
    Fs = round(1/(timeVals(2)-timeVals(1)));
    freq = 0:1/diff(timeRange):Fs-1/diff(timeRange);
    % finding positions of data points that fall between the selected
    % timeRange (stimulus and baseline are two different stimulus types therefor no stPos and blPos separately)
    Pos = timeVals >= timeRange(1) & timeVals < timeRange(2); 

    %% doing fourier transform and calculating power
   
    for icond = 1:length(types)
        if icond>1
            if useERP == 1
                fftStim{icond}{ielec,1} = (abs(fft(mean(analogDataAllGra{icond}{ielec,1}(:,Pos),1)))/length(freq));
                powStim{icond}{ielec,1} = fftStim{icond}{ielec,1}.^2;
            else
                fftStim{icond}{ielec,1} = (mean(abs(fft(analogDataAllGra{icond}{ielec,1}(:,Pos),[],2)))/length(freq));
                powStim{icond}{ielec,1} = fftStim{icond}{ielec,1}.^2;
            end
        else 
            for rfi = 1:length(unique(RfVals))
                for rvi = 1:length(unique(RvVals))
                    if useERP == 1
                        fftStim{icond}{ielec,rfi,rvi} = (abs(fft(mean(analogDataAllGra{icond}{ielec,rfi,rvi}(:,Pos),1)))/length(freq));
                        powStim{icond}{ielec,rfi,rvi} = fftStim{icond}{ielec,rfi,rvi}.^2;
                    else
                        fftStim{icond}{ielec,rfi,rvi} = (mean(abs(fft(analogDataAllGra{icond}{ielec,rfi,rvi}(:,Pos),[],2)))/length(freq));
                        powStim{icond}{ielec,rfi,rvi} = fftStim{icond}{ielec,rfi,rvi}.^2;
                    end
                end
            end
            
        end
    end
end

TrialSortedData= analogDataAllGra;
TrialAverageData = avgAnalogData; %%% data averaged across trials

%%%%%% averaging power across electrodes %%%%%%
if  strcmpi(referenceType,'bipolarRef')
    goodElecs = ElecSet;
else
    goodElecs = cell2mat(good_elecs{index});
end


for icond = 1:length(powStim)
    if icond>1
        avgPow_acrssElecs{icond} = mean(cell2mat(powStim{1,icond}(goodElecs)));
    else
        for rf = 1:length(uRf)
            for rv = 1:length(uRv)
                avgPow_acrssElecs{icond}{rf,rv}= mean(cell2mat(powStim{1,icond}(goodElecs,rf,rv)));
            end
        end
    end
end

%% time frequency analysis
movingwin = [0.3 0.025];
%%% params %%%
if ~exist('params','var') || isempty(params)
    params.tapers   = [1 1];
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 100];
    params.trialave = 1; %averaging across trials
end

clear Sbl Sbl_meanElec  Sst Sst_meanElec SstR SstR_meanElec
for iele = 1:length(goodElecs)
    for iRf = 1:length(uRf)
        for iRv = 1:length(uRv)
            [SstR{iRf,iRv}(iele,:,:,:),T,FreqVals]=mtspecgramc(TrialSortedData{1,1}{goodElecs(iele),iRf,iRv}',movingwin,params);
        end
    end
    [Sst(iele,:,:,:),~,~]=mtspecgramc(TrialSortedData{1,2}{goodElecs(iele),1}',movingwin,params);
    [Sbl(iele,:,:,:),~,~]=mtspecgramc(TrialSortedData{1,3}{goodElecs(iele),1}',movingwin,params);
end

TimeVals = T+timeVals(1)-1/Fs;
%averging across elecs
Sbl_meanElec = squeeze(mean(Sbl,1)); % averaging across elecs - baseline condition
Sst_meanElec = squeeze(mean(Sst,1)); %  averaging across elecs - amplitude modulated sinusoid condition
SstR_meanElec =cellfun(@squeeze,(cellfun(@mean,SstR,'UniformOutput',false)),'UniformOutput',false); %ripple conditions
Sbl_meanElec_rep = repmat({Sbl_meanElec},length(uRf),length(uRv));
powChange{1,1} = cellfun(@rdivide,SstR_meanElec,Sbl_meanElec_rep,'Un',0); % dividing 
powChange{1,2} = Sst_meanElec./Sbl_meanElec;

end