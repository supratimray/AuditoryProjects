%%%%% stimulus set for experiment one aim2 %%%%%%%

clear
close all
clc

modIndex = 0.9;
modFreq = 40;

count = 1;
for mf = 1:length(modFreq)
    %figure,
    for m = 1:length(modIndex)
        
        fs =44100; % sampling frequency
        t=0:1/fs:0.8; % Total time for simulation
        
        %XXXXXXXXXXXXXXXXXXXXX carrier signal generation XXXXXXXXXXXXXXXXXXXXXXXXXX
        Ac=1;% Amplitude of carrier signal 
        fc=1000;% Frequency of carrier signal
        Tc=1/fc;% Time period of carrier signal
        yc=Ac*sin(2*pi*fc*t);% Equation of carrier signal
        
        Am=Ac*modIndex(m);%[ where, modulation Index (m)=Am/Ac ] % Amplitude of modulating signal
        fa=modFreq(mf); % Frequency of modulating signal
        Ta=1/fa;% Time period of modulating signal
        ym=Am*sin(2*pi*fa*t); % Equation of modulating signal
        
        
        
        %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX AM Modulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        y_notnorm=Ac*(1+modIndex(m)*sin(2*pi*fa*t)).*sin(2*pi*fc*t);% Equation of Amplitude
        %modulated signal
        y = y_notnorm./max(abs(y_notnorm));
        mindex = modIndex(m);
        yall(count,:) = y;
%         subplot(3,2,count);
%         plot(t,y), grid on; %ylim([-10 10])% Graphical representation of Modulating signal
%         title ( '  Modulating Signal   ');
%         xlabel ( ' time(sec) ');
%         ylabel (' Amplitude(volt)   ');
         count = count +1;
        str = sprintf('Azi_0.0_Elev_0.0_Type_2_CF_1000_MF_40_MD_0.9_Dur_800.wav');
        audiowrite(str,y,fs)
        
    end
end
