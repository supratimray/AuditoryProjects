%%%%%%%%%% generating conditions file for for experiment 1  %%%%%%%%%%%%
clear; close all;clc;
%%%% opening conditions file %%%%

fid =fopen('Aud_Gamma_conditions.txt','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t','Condition','Frequency','Block','Timing File','TaskObject#1');

% for generating sound file name

%%ripple stimuli
azi = 0;
ele = 0;
type = 1;
rippleFreq = [0,0.4,0.8,1.6,2.4];
rp = 0;
md = 0.9;
rippleVel = [0,2.5,5,10,20,40];

count = 1;
for rf = 1:length(rippleFreq)
    for rv = 1:length(rippleVel)
        str = ['snd(Azi_' sprintf('%.1f',azi) '_Elev_' sprintf('%.1f',ele) '_Type_' num2str(type) '_RF_' sprintf('%.1f',rippleFreq(rf)) '_RP_' num2str(rp) '_MD_' sprintf('%.1f',md) '_RV_' sprintf('%.1f',rippleVel(rv)) '_Dur_' num2str(800) '.wav)'];
        fprintf(fid,'\n%d\t%d\t%d\t%s\t%s',count,1,1,'timingfile_humanECoG', str);
        count = count+1;
    end
end
 
%%ASSR 
str = sprintf('snd(Azi_0.0_Elev_0.0_Type_2_CF_1000_MF_40_MD_0.9_Dur_800.wav)');
fprintf(fid,'\n%d\t%d\t%d\t%s\t%s',count,1,1,'timingfile_humanECoG', str);

%%blank stimuli
fprintf(fid,'\n%d\t%d\t%d\t%s\t%s\t%s\t',count+1,1,1,'timingfile_humanECoG','snd(Noise_Dur_800.wav)');
fclose(fid);

