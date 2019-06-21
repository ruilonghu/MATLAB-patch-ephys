%% cell attached

tic
close all
clearvars -except Trace_*

%% settings and shit
savename = '190619 slice1 (1) AOB GC';
yes_save = 1;

% which experiment
Experimentnum = 1;
% load which trace?
Trace = 1;
channel = 1;
channel_B = [];


% threshold after the high pass filter
threshold = (10^-12) * 4; % 


% yaxis off?
yOFF = 1;
ypA_scalebar = 10; % pA for cale bar


% drug info
drugpresent = 0;
drugname = 'Bic 10uM';
drugON = 60 * 10; % seconds (from minutes)
drugOFF = 60 * 15; % seconds

% do autocorrelation?
autoCORR = 1;

%%  set up the list of traces and find Fs
Experiment = ['Trace_' num2str(Experimentnum) '_'];
tracelist = who([Experiment, num2str(Trace), '_*_', num2str(channel)]);

for i = 1:length(tracelist)
    reord_tracelist_I{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(channel)];
end

if isempty(channel_B) ~= 1;
    for i = 1:length(tracelist)
        reord_tracelist_I_b{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(channel_B)];
    end
end

tmp = eval(reord_tracelist_I{1});
maxT = tmp(end,1); % time index for end
Fs = length(tmp(:,1))/maxT;


%% stack into long trace and do stuff

longT = [];
longI = [];
lastTj = 0;

for j = 1:length(reord_tracelist_I)
    tmp = eval(reord_tracelist_I{j});
    tmpdata_T = tmp(:,1);
    tmpdata_V = tmp(:,2);
    longT = vertcat(longT, (tmpdata_T + lastTj));
    longI = vertcat(longI, tmpdata_V) ;
    lastTj = longT(end); % it'll be the end time for hte next loop iteration
    
end


if isempty(channel_B) ~= 1;
    
    longT_B = [];
    longI_B = [];
    lastTj = 0;
    for j = 1:length(reord_tracelist_I_b)
        tmp = eval(reord_tracelist_I_b{j});
        tmpdata_T = tmp(:,1);
        tmpdata_V = tmp(:,2);
        longT_B = vertcat(longT_B, (tmpdata_T + lastTj));
        longI_B = vertcat(longI_B, tmpdata_V) ;
        lastTj = longT(end); % it'll be the end time for hte next loop iteration
        
    end
end

%% high pass filter and detect spikes
[bLP,aLP] = butter(1, 500/(Fs/1),'high'); %300Hz LP filter 2nd order
% default butter(2, 200/(Fs/2), 'low');

longIf = (filtfilt(bLP,aLP,longI)); % Low pass @100Hz FOR SUPER FILTER
longIf_inv = -(filtfilt(bLP,aLP,longI));

if isempty(channel_B) ~= 1;
    longI_Bf = (filtfilt(bLP,aLP,longI_B));
    longI_Bf_inv = -(filtfilt(bLP,aLP,longI_B));% Low pass @100Hz FOR SUPER FILTER
end

% detect spikes

[pks, lct] = findpeaks(longIf_inv, 'MinPeakHeight',threshold, 'MinPeakDistance', round(.001*Fs));
PSC_eventT = longT(lct);

if isempty(channel_B) ~= 1;
    [pks, lct] = findpeaks(longI_Bf_inv, 'MinPeakHeight',threshold, 'MinPeakDistance', round(.001*Fs));
    PSC_eventT_B = longT_B(lct);
end

%% do autocorrelation
if autoCORR == 1
    combotemp=PSC_eventT;
    maxtime=max(combotemp);
    maxtime=ceil(maxtime);%Times are in seconds
    timebins=0:1:maxtime;
    TbinsA=histc(PSC_eventT,timebins);
    % TbinsB=histc(PSC_eventT_B,timebins);
    maxlag= 60;
    [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
    
    
    figure; plot(plotlags,c1, 'k');
    xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    
    box off;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','inches');
    set(gcf,'Position',[1 1 10 7.5]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 7.5]);
    set(gcf, 'PaperPosition', [0 0 10 7.5]);
    set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
    
    norm_histogramA = gca;
    
    if isempty(channel_B) ~= 1;
        combotemp=PSC_eventT_B;
        maxtime=max(combotemp);
        maxtime=ceil(maxtime);%Times are in seconds
        timebins=0:1:maxtime;
        TbinsA=histc(PSC_eventT_B,timebins);
        % TbinsB=histc(PSC_eventT_B,timebins);
        maxlag=60;
        [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
        
        
        figure; plot(plotlags,c1, 'k');
        xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        ylabel('Normalized Event Auto-Correlation cell 2', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        
        box off;
        set(gcf,'Color',[1 1 1]);
        set(gcf,'Units','inches');
        set(gcf,'Position',[1 1 10 7.5]);
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [10 7.5]);
        set(gcf, 'PaperPosition', [0 0 10 7.5]);
        set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
        
        norm_histogramB = gca;
    end
    
    
    
    % do cross correlation
    if isempty(channel_B) ~= 1;
        combotemp=PSC_eventT;
        maxtime=max(combotemp);
        maxtime=ceil(maxtime); %Times are in seconds
        timebins=0:.1:maxtime;
        TbinsA=histc(PSC_eventT,timebins);
        TbinsB=histc(PSC_eventT_B,timebins);
        maxlag=60;
        [c1, plotlags] = xcorr(TbinsA,TbinsB,maxlag,'coeff');
        
        figure; plot(plotlags,c1, 'k');
        xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        ylabel('Normalized Event Cross-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        
        box off;
        set(gcf,'Color',[1 1 1]);
        set(gcf,'Units','inches');
        set(gcf,'Position',[1 1 10 7.5]);
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [10 7.5]);
        set(gcf, 'PaperPosition', [0 0 10 7.5]);
        set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
        
        xcorr_hist = gca;
    end
    
end


%% calculate spike rate
binT = 1; % .1 s bin size
binIwin = [1:binT:longT(end)];

for i = 1:length(binIwin)
    if i == length(binIwin)
        Tstop = longT(end);
    else
        Tstop = binIwin(i+1);
    end
    Tstart = binIwin(i);
    spk_bin_ct(i) = sum(length(find(PSC_eventT >  Tstart & PSC_eventT < Tstop)));
    spk_bin_r(i) = spk_bin_ct(i)./binT;
end

if isempty(channel_B) ~= 1;
    binT = 1; % .1 s bin size
    binIwin_B = [1:binT:longT_B(end)];
    
    for i = 1:length(binIwin)
        if i == length(binIwin)
            Tstop = longT_B(end);
        else
            Tstop = binIwin(i+1);
        end
        Tstart = binIwin(i);
        spk_bin_ct_B(i) = sum(length(find(PSC_eventT_B >  Tstart & PSC_eventT_B < Tstop)));
        spk_bin_r_B(i) = spk_bin_ct_B(i)./binT;
    end
end

%% plots IPSCs
figure; subplot (3,1,[1 2]);
plot(longT, longIf, 'k')
hold on


    yvect = [(mean(longI) + 8*std(longI)) (mean(longI) + 6*std(longI))];

for i = 1:length(PSC_eventT)
    plot([PSC_eventT(i), PSC_eventT(i)], yvect, 'r')
end

if drugpresent == 1
    drugvect = [drugON drugOFF];
    

yvaldrug = mean(longI) + 15*std(longI);

    
    yvect = [yvaldrug yvaldrug];
    plot(drugvect, yvect, 'b', 'LineWidth', 3)
    text(mean(drugvect), (yvaldrug + std(longI)), drugname, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
    % calculate drug times too
    dataPSC.pre_drugR = length(find(PSC_eventT < drugON))/(drugON); % PSC per second
    dataPSC.on_drugR = length(find(PSC_eventT > drugON & PSC_eventT < drugOFF))/(drugOFF - drugON); % PSC per second
    
end

box off;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 10 7.5]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 7.5]);
set(gcf, 'PaperPosition', [0 0 10 7.5]);
if yOFF == 1
    set(gca, 'YTick', []);
    set(gca,'YColor','w');
end
set(gca, 'XTick', 0:60:longT(end), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
xlabel('Time (s)');
xlim([0 longT(end)])
% general axis limits

    ylim([(mean(longI) -6*std(longI)) (mean(longI) + 10*std(longI))])

% plot pA scale bar
yscale = (10^-12) * ypA_scalebar;


    yvect = [(mean(longI) - 2*std(longI)) (mean(longI) - 2*std(longI) + yscale)];

    
plot([1, 1],yvect, 'k', 'linewidth',3)
text(5, mean(yvect), [num2str(ypA_scalebar) ' pA' ], 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')

title([savename], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');


%plots firing rate
subplot (3,1,3);
plot(binIwin, spk_bin_r, 'k')
set(gca, 'XTick', 0:60:longT(end), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
xlabel('Time (s)');
ylabel('Instantaneous Firing Rate (Hz)','Fontsize',10);
xlim([0 longT(end)])

IPSCrasterA = gca;



%%
if isempty(channel_B) ~= 1;
    figure; subplot (3,1,[1 2]);
    plot(longT_B, longI_B, 'k')
    hold on
    if EPSCs == 1
        yvect = [(mean(longI_B) + 10*std(longI_B)) (mean(longI_B) + 8*std(longI_B))];
    else
        yvect = [(mean(longI_B) - 11*std(longI_B)) (mean(longI_B) - 9*std(longI_B))];
    end
    for i = 1:length(PSC_event_B)
        plot([PSC_eventT_B(i), PSC_eventT_B(i)], yvect, 'r')
    end
    
    if drugpresent == 1
        drugvect = [drugON drugOFF];
        
        if EPSCs == 1
            yvaldrug = mean(longI_B) + 15*std(longI_B);
        else
            yvaldrug = mean(longI_B) + 28*std(longI_B);
        end
        
        yvect = [yvaldrug yvaldrug];
        plot(drugvect, yvect, 'b', 'LineWidth', 3)
        text(mean(drugvect), (yvaldrug + std(longI_B)), drugname, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
        % calculate drug times too
        dataPSC.pre_drugR = length(find(PSC_eventT_B < drugON))/(drugON); % PSC per second
        dataPSC.on_drugR = length(find(PSC_eventT_B > drugON & PSC_eventT_B < drugOFF))/(drugOFF - drugON); % PSC per second
        
    end
    
    box off;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','inches');
    set(gcf,'Position',[1 1 10 7.5]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 7.5]);
    set(gcf, 'PaperPosition', [0 0 10 7.5]);
    if yOFF == 1
        set(gca, 'YTick', []);
        set(gca,'YColor','w');
    end
    set(gca, 'XTick', 0:60:longT_B(end), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    xlabel('Time (s)');
    xlim([0 longT_B(end)])
    % general axis limits


        ylim([(mean(longI_B) -13*std(longI_B)) (mean(longI_B) + 30*std(longI_B))])
  
    
    % plot pA scale bar
    yscale = (10^-12) * ypA_scalebar;
    
    if EPSCs == 1
        yvect = [(mean(longI_B) - 25*std(longI_B)) (mean(longI_B) - 25*std(longI_B) + yscale)];
    else
        yvect = [(mean(longI_B) - 8*std(longI_B)) (mean(longI_B) - 8*std(longI_B) + yscale)];
    end
    plot([1, 1],yvect, 'k', 'linewidth',3)
    text(5, mean(yvect), [num2str(ypA_scalebar) ' pA' ], 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
    
    title([savename], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
    
    %plots firing rate
    subplot (3,1,3);
    plot(binIwin_B, spk_bin_r_B, 'k')
    set(gca, 'XTick', 0:60:longT_B(end), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Instantaneous Firing Rate (Hz)','Fontsize',10);
    xlim([0 longT_B(end)])
    
    IPSCrasterB = gca;
end



if yes_save == 1
    saveas(IPSCrasterA, [savename ' cell-attached raster.jpg'], 'jpg')
    if autoCORR ==1
        saveas(norm_histogramA, [savename ' cell-attached normalized autocorrelation.jpg'], 'jpg')
    end
    
    if isempty(channel_B) ~= 1;
        saveas(IPSCrasterB, [savename ' 2 cell-attached raster.jpg'], 'jpg')
        saveas(xcorr_hist, [savename ' cell-attached normalized cross-correlation.jpg'], 'jpg')
        if autoCORR ==1
            saveas(norm_histogramB, [savename ' cell-attached normalized autocorrelation.jpg'], 'jpg')
        end
    end
    
    save([savename, ' cellattached_FR.mat'], '-regexp',  '^(?!Trace_.*$).')

end

total_time = toc