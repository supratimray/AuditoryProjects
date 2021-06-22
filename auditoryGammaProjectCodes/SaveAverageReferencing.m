%%%%%% creating and saving vals for average referencing %%%%%
clear;clc;
[subjectName,expDate,protocolNames,good_elecs] = EcogAuditoryGammaData;
folderSourceString = 'D:\OneDrive - Indian Institute of Science\divya\NimhansRippleProject\Divya_AuditoryProjects\data\humanECoG';
gridType = 'ECoG';
%%
for id = 1:length(subjectName)
    
    % Get folders
    folderName = string(fullfile(folderSourceString,'data',subjectName{id},gridType,expDate{id},protocolNames{id}));
    folderExtract = fullfile(folderName,'extractedData');
    folderSegment = fullfile(folderName,'segmentedData');
    folderLFP = fullfile(folderSegment,'LFP');
   
   
    x = load(fullfile(folderLFP,'lfpInfo.mat'));
    ElecSet=sort(x.analogChannelsStored);
    
    
    all_elec_Data = [];
    for i = 1:length(ElecSet)-1%elecs
        all_Data = load(fullfile(folderLFP, ['elec' num2str(ElecSet(i)) '.mat']));
        all_elec_Data = [all_elec_Data;all_Data.analogData];
    end
    
    analogData= mean(all_elec_Data);
    analogData= repmat(analogData,320,1);
    save(fullfile(folderLFP,'\AvgRef.mat'),'analogData');
end