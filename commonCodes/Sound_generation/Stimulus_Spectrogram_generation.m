
%%%%%%% to create spectrogram of stimuli set %%%%%
RF = [0.0,0.4,0.8,1.6,2.4];

RV = [0,2.5,5,10,20,40];
for i = 0:3
    X = 2^i*5;
    RV = [RV,X];
end
data = cell(1,32);

count = 1;
for iRF = 1:length(RF)
    for iRV = 1:length(RV)
        rf = sprintf('%0.1f',RF(iRF));     
        rv = sprintf('%0.1f',RV(iRV));
        A = ("Azi_0.0_Elev_0.0_Type_1_RF_"+rf+"_RP_0_MD_0.9_RV_" + rv + "_Dur_800.wav");
        [data{1,count},fs{count}] = audioread (A);
        count = count+1;
    end
end

[data{1,31},fs{1,31}] = audioread('Azi_0.0_Elev_0.0_Type_2_CF_1000_MF_40_MD_0.9_Dur_800.wav');
[data{1,32},fs{1,32}] = audioread('Noise_Dur_800.wav');

figure,
[hplot,~] = getPlotHandles(6,5,[0.05 0.09 0.65 0.9],0.01,0.005,1);
h2 = hplot';
for ispec = 1:30
    y = data{1,ispec}(:,1);
    Nx = length(y);
    nsc = floor(Nx/100);
    nov = floor(nsc/3);
    nff = 8192;
    subplot(h2(ispec))
    spectrogram(y,hanning(nsc),nov,nff,44100,'yaxis');
    ylim([0 8]);
    colormap('jet');
    colorbar('off');
    set(gca,'xlabel',[]);
    set(gca,'ylabel',[]);
    set(gca,'fontsize',12)
%     if ispec == 21
%         xlabel('Time (ms)');
%         ylabel('Frequency (kHz)');
%     end
%     if ispec>21 || ispec<21
%         set(gca,'xticklabels',[]);
%         set(gca,'yticklabels',[]);
%     end
end

subplot('Position',[0.73 0.6 0.25 0.3])
spectrogram(data{1,31}(:,1),hanning(nsc),nov,nff,44100,'yaxis');
ylim([0 2]);
colormap('jet');
colorbar('off');
set(gca,'fontsize',11);

subplot('Position',[0.73 0.2 0.25 0.3])
spectrogram(data{1,32}(:,1),hanning(nsc),nov,nff,44100,'yaxis');
set(gca,'fontsize',11);colorbar('off');