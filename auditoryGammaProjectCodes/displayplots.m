%This code generates PSD and change in power spectrogram for all different
%conditions 
% Input includes - subjectId, 
% AvgPow_acrssElecs - Array of power during different stimulus
% conditions which is averaged across electrodes. 
%powChange - Array of power as a function of both time and frequency which
%is averaged across electrodes
% TimeVals - Vector indicating timepoints at which sample was collected
% FreqVals - Frequency vector
%uRf - Unique Ripple frequencies
%uRv - Unique ripple velocity
%timeRange = time period for which PSD is to be plotted

function displayplots(subjectIndex,AvgPow_acrssElecs,powChange,TimeVals,FreqVals,uRf,uRv,timeRange,Samp_Freq)

if ~exist('timeRange','var');timeRange = [0.25 0.75];end
if ~exist('Samp_Freq','var');Samp_Freq = 2500; end
% plotting power averaged across electrodes %

% making frequency axis
Fs = Samp_Freq; %sampling frequency
Freq = 0:1/diff(timeRange):Fs-1/diff(timeRange);
%%
figure,
[hplot,~] = getPlotHandles(length(uRf),length(uRv),[0.04 0.065 0.65 0.84],0.005,0.035,1);
h2 = hplot';
xtcks = (0:30:90);
count = 1;
for subjectid = subjectIndex
    for icond = 1
        for jrf = 1:length(uRf)
            for jrv = 1:length(uRv)
                subplot(h2(count))
                plot(Freq,10*log10(AvgPow_acrssElecs{1,subjectid}{1,icond}{jrf,jrv}),'r','Linewidth',2.5); hold on; plot(Freq,10*log10(AvgPow_acrssElecs{1,subjectid}{1,3}),'g','Linewidth',2.5)
                xlim([0 100]);title(uRf(jrf) + "cyc/oct  " + uRv(jrv) + "cyc/sec" )
                count = count+1;
                if jrv == 1
                    ylabel('Power (dB)');
                else
                    set(gca,'yticklabels',[]);
                end
                if count >25
                    xlabel('Frequency (Hz)');
                    xticks(xtcks);
                else
                    set(gca,'xticklabels',[]);
                end
            end
            legend('Stim','Base');
        end
    end
end

sgtitle("Power Spectral Densities of Stimulus and baseline conditions averaged across electrodes   subjectId" + subjectid')
subplot('Position',[0.73 0.6 0.25 0.3])
plot(Freq,10*log10(AvgPow_acrssElecs{1,subjectid}{1,2}),'r','Linewidth',2.5); hold on; plot(Freq,10*log10(AvgPow_acrssElecs{1,subjectid}{1,3}),'g','Linewidth',2.5)
 xlim([0 100]);title('Amplitude Modulated Sinusoid at 40Hz');
 ylabel('Power (dB)');legend('Stim','Base');xlabel('Frequency (Hz)');
 
 


%%%%%%%% change in power %%%%%%%%%%

figure,
[hplot,~] = getPlotHandles(length(uRf),length(uRv),[0.04 0.065 0.65 0.84],0.005,0.035,1);
h2 = hplot';

count = 1;
for subjectid = subjectIndex
    for icond = 1
        for jrf = 1:length(uRf)
            for jrv = 1:length(uRv)
                subplot(h2(count))
                pcolor(TimeVals{subjectid}, FreqVals{subjectid}, 10*log10(powChange{1,subjectid}{1,icond}{jrf,jrv}'));
                colormap('jet'); caxis([-10 10]); 
                shading interp;
                title(uRf(jrf) + "cyc/oct  " + uRv(jrv) + "cyc/sec" )
                count = count+1;
                if jrv == 1
                    ylabel('Frequency (Hz)','FontWeight','bold');
                else
                    set(gca,'yticklabels',[]);
                end
                if jrv == 6
                    a=colorbar; 
                    a.Label.String = 'Power (dB)';
                    a.FontWeight = 'bold';
                    a.Position = [0.7 0.065 0.01 0.149 ];
                end
                if count >25
                    xlabel('Time(sec)','FontWeight','bold');
                else
                    set(gca,'xticklabels',[]);
                end
            end
        end
    end
end
sgtitle("Change in Power during stimulus and baseline conditions (averaged across electrodes)    subjectId" + subjectid)

subplot('Position',[0.73 0.4 0.25 0.3])
pcolor(TimeVals{subjectid}, FreqVals{subjectid}, 10*log10(powChange{1,subjectid}{1,icond+1}'));
colormap('jet');set(gca,'FontWeight','bold')
caxis([-10 10]);
shading interp;
a=colorbar;
a.Label.String = 'Power (dB)';
xlabel('Time(sec)');ylabel('Frequency (Hz)');
title('Amplitude modulated sinusoid at 40Hz');

end