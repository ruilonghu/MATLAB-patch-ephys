tic
close all
clearvars -except Trace_*

%% settings and shit
savename = '190409 slice1 (1) MOB dSAC - NA';
yes_save = 1;
% which experiment
Experimentnum =        1        ;
% load which trace
Trace =                9        ;
channel =              1        ;
channel_B = [];

% Threshold for spikes
threshold = -.0150;

% yaxis off?
yOFF = 1;
% ymV_scalebar = 10; % mV for scale bar

scaleX1 = 2 * 60; % 2 min
scaleY = 20 * 10^-3; % 10 mV

% drug text offset
txtoffset_A = .05; % default .1
txtoffset_B = .1;
txtoffset_C = .2;


% drug info
drugs.drugpresent = 1;
drugs.drugname = 'NA 10uM'; %subP 500nM
drugs.drugON = 60 *  61; % seconds (from minutes)
drugs.drugOFF = 60 * 63; % seconds

get_number_drugon_only = 1; % do you want to grab numbers from only when the drug was on (1), basically adding 1 minute to drugON and drugOFF.
% otherwise, it do everything until the trace ends (from when drug was turned on)
drugs.drugonDELAY = 60*3;
drugs.drugonANAL = 60*3;

drugs.drug2present = 1;
drugs.drugname2 = 'muscarine 10uM';
drugs.drugON2 = 60 *  45; % seconds (from minutes)
drugs.drugOFF2 = 60 * 50; % seconds

drugs.drug3present =1;
drugs.drugname3 = '4DAMP 100 nM';
drugs.drugON3 = 60 * 38; % seconds (from minutes)
drugs.drugOFF3 = 60 * 55 % seconds


% do autocorrelation?
autoCORR = 0;
msbin = 0;

%%  set up the list of traces and find Fs
Experiment = ['Trace_' num2str(Experimentnum) '_'];
tracelist = who([Experiment, num2str(Trace), '_*_', num2str(channel)]);

for i = 1:length(tracelist)
    reord_tracelist_V{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(channel)];
end

if isempty(channel_B) ~= 1
    for i = 1:length(tracelist)
        reord_tracelist_V_b{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(channel_B)];
    end
end

tmp = eval(reord_tracelist_V{1});
maxT = tmp(end,1); % time index for end
Fs = length(tmp(:,1))/maxT;


%% stack into long trace and do stuff

for i = 1:length(reord_tracelist_V)
    tmp = eval(reord_tracelist_V{i});
    maxT = tmp(end,1);
    maxidx = length(tmp(:,1));
    idxstart = ((i-1)*maxidx) + 1;
    idxstop = (i*maxidx);
    longT(idxstart:idxstop) = (tmp(:,1) + ((i-1)*maxT));
    longV(idxstart:idxstop) = tmp(:,2);
end

if isempty(channel_B) ~= 1
    for i = 1:length(reord_tracelist_V_b)
        tmp = eval(reord_tracelist_V_b{i});
        maxT = tmp(end,1);
        maxidx = length(tmp(:,1));
        idxstart = ((i-1)*maxidx) + 1;
        idxstop = (i*maxidx);
        longT_B(idxstart:idxstop) = (tmp(:,1) + ((i-1)*maxT));
        longV_B(idxstart:idxstop) = tmp(:,2);
    end
end

%% detect spikes

[pks, lct] = findpeaks(longV, 'MinPeakHeight',threshold, 'MinPeakDistance', round(.001*Fs));
spk_eventT = longT(lct);

%%
if isempty(channel_B) ~= 1;
    [pks, lct] = findpeaks(longV_B, 'MinPeakHeight',threshold, 'MinPeakDistance', round(.001*Fs));
    spk_eventT_B = longT_B(lct);
end

%% calculate spike rate
binT = 1; % 1 s bin size
binIwin = [1:binT:longT(end)];

for i = 1:length(binIwin)
    if i == length(binIwin)
        Tstop = longT(end);
    else
        Tstop = binIwin(i+1);
    end
    Tstart = binIwin(i);
    spk_bin_ct(i) = length(find(spk_eventT >  Tstart & spk_eventT < Tstop));
    spk_bin_r(i) = spk_bin_ct(i)./binT;
end


if isempty(channel_B) ~= 1;
    binIwinB = [1:binT:longT_B(end)];
    for i = 1:length(binIwinB)
        if i == length(binIwinB)
            Tstop = longT_B(end);
        else
            Tstop = binIwinB(i+1);
        end
        Tstart = binIwinB(i);
        spk_bin_ct_B(i) = length(find(spk_eventT_B >  Tstart & spk_eventT_B < Tstop));
        
        spk_bin_r_B(i) = spk_bin_ct_B(i)./binT;
        
        
    end
end


%% do autocorrelation
% if autoCORR == 1
% combotemp=spk_eventT;
% maxtime=max(combotemp);
% maxtime=ceil(maxtime);%Times are in seconds
% timebins=0:.001:maxtime;
% TbinsA=histc(spk_eventT,timebins);
% % TbinsB=histc(PSC_eventT_B,timebins);
% maxlag= 500;
% [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
%
%
% figure; plot(plotlags,c1, 'k');
% xlabel('Time Lag in 1ms bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
% ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
%
% box off;
% set(gcf,'Color',[1 1 1]);
% set(gcf,'Units','inches');
% set(gcf,'Position',[1 1 10 7.5]);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [10 7.5]);
% set(gcf, 'PaperPosition', [0 0 10 7.5]);
% set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
%
% norm_histogram = gca;
%
% if isempty(channel_B) ~= 1;
% combotemp=spk_eventT_B;
% maxtime=max(combotemp);
% maxtime=ceil(maxtime);%Times are in seconds
% timebins=0:1:maxtime;
% TbinsA=histc(spk_eventT_B,timebins);
% % TbinsB=histc(PSC_eventT_B,timebins);
% maxlag=60;
% [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
%
%
% figure; plot(plotlags,c1, 'k');
% xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
% ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
%
% box off;
% set(gcf,'Color',[1 1 1]);
% set(gcf,'Units','inches');
% set(gcf,'Position',[1 1 10 7.5]);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [10 7.5]);
% set(gcf, 'PaperPosition', [0 0 10 7.5]);
% set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
%
% norm_histogram_B = gca;
%
% end
%
%
% end

if autoCORR == 1
    
    if drugs.drugpresent == 0
        combotemp=spk_eventT;
        maxtime=max(combotemp);
        maxtime=ceil(maxtime);%Times are in seconds
        
        
        if msbin == 1
            tbinwin = .001;
            maxlag = 1000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins=0:tbinwin:maxtime;
        TbinsA=histc(spk_eventT,timebins);
        %     maxlag=20;
        [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
        
        
        figure; plot(plotlags,c1, 'k');
        
        if msbin == 1
            xlabel('Time Lag in 1ms bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        else
            xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        end
        
        ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
        
        box off;
        set(gcf,'Color',[1 1 1]);
        set(gcf,'Units','inches');
        set(gcf,'Position',[1 1 10 7.5]);
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [10 7.5]);
        set(gcf, 'PaperPosition', [0 0 10 7.5]);
        set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
        
        norm_histogram = gca;
        
        if yes_save == 1
            saveas(norm_histogram, [savename ' normalized autocorrelation.jpg'], 'jpg')
        end
        
        
        % if drug is present
    elseif drugs.drugpresent == 1
        
        % pre drug
        combotemp = spk_eventT(1:length(find(spk_eventT < drugs.drugON)));
        
        maxtime=max(combotemp);
        maxtime=ceil(maxtime);%Times are in seconds
        
        
        if msbin == 1
            tbinwin = .001;
            maxlag = 1000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins=0:tbinwin:maxtime;
        TbinsA=histc(spk_eventT,timebins);
        %     maxlag=20;
        [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
        
        if ~isempty(combotemp)
            figure; plot(plotlags,c1, 'k');
            
            if msbin == 1
                xlabel('Time Lag in 1ms bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
            else
                xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
            end
            
            ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
            title(['PRE DRUG'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
            
            box off;
            set(gcf,'Color',[1 1 1]);
            set(gcf,'Units','inches');
            set(gcf,'Position',[1 1 10 7.5]);
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [10 7.5]);
            set(gcf, 'PaperPosition', [0 0 10 7.5]);
            set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
            
            norm_histogram_pre = gca;
            
            
            if yes_save == 1
                saveas(norm_histogram_pre, [savename ' normalized autocorrelation pre-drug.jpg'], 'jpg')
            end
        end
        
        
        
        % while drug is on
        
        if get_number_drugon_only == 1
            drugdrugTMP = find(spk_eventT > (drugs.drugON + drugs.drugonDELAY) & spk_eventT < (drugs.drugOFF + drugs.drugonDELAY));
        elseif get_number_drugon_only == 0
            drugdrugTMP = find(spk_eventT > (drugs.drugON + drugs.drugonDELAY));
        end
        
        combotemp = drugdrugTMP;
        
        maxtime=max(combotemp);
        maxtime=ceil(maxtime);%Times are in seconds
        
        if msbin == 1
            tbinwin = .001;
            maxlag = 1000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins=0:tbinwin:maxtime;
        TbinsA=histc(spk_eventT,timebins);
        %     maxlag=20;
        [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
        
        if ~isempty(combotemp)
            
            figure; plot(plotlags,c1, 'k');
            
            if msbin == 1
                xlabel('Time Lag in 1ms bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
            else
                xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
            end
            
            ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
            title(['ON DRUG'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
            
            box off;
            set(gcf,'Color',[1 1 1]);
            set(gcf,'Units','inches');
            set(gcf,'Position',[1 1 10 7.5]);
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperSize', [10 7.5]);
            set(gcf, 'PaperPosition', [0 0 10 7.5]);
            set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
            
            norm_histogram_post = gca;
            
            
            if yes_save == 1
                saveas(norm_histogram_post, [savename ' normalized autocorrelation ON-drug.jpg'], 'jpg')
            end
        end
        
        % drug washout
        
        if get_number_drugon_only == 1
            drugdrugTMP = find(spk_eventT > (longT(end) - drugs.drugonANAL) & spk_eventT < (longT(end)));
            
            combotemp = drugdrugTMP;
            
            maxtime=max(combotemp);
            maxtime=ceil(maxtime);%Times are in seconds
            
            if msbin == 1
                tbinwin = .001;
                maxlag = 1000; %ms
            else
                tbinwin = 1;
            maxlag = 60; %s
            end
            
            timebins=0:tbinwin:maxtime;
            TbinsA=histc(spk_eventT,timebins);
            %     maxlag=20;
            [c1, plotlags] = xcorr(TbinsA,maxlag,'coeff');
            
            if ~isempty(combotemp)
                
                figure; plot(plotlags,c1, 'k');
                
                if msbin == 1
                    xlabel('Time Lag in 1ms bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
                else
                    xlabel('Time Lag in 1s bin', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
                end
                
                ylabel('Normalized Event Auto-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
                title(['Washout DRUG'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
                
                box off;
                set(gcf,'Color',[1 1 1]);
                set(gcf,'Units','inches');
                set(gcf,'Position',[1 1 10 7.5]);
                set(gcf, 'PaperUnits', 'inches');
                set(gcf, 'PaperSize', [10 7.5]);
                set(gcf, 'PaperPosition', [0 0 10 7.5]);
                set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
                
                norm_histogram_washout = gca;
                
                if yes_save == 1
                    saveas(norm_histogram_washout, [savename ' normalized autocorrelation washout-drug.jpg'], 'jpg')
                end
            end
        end
        
        
    end
end
%% Do cross correlation
if isempty(channel_B) ~= 1;
    combotemp=horzcat(spk_eventT,spk_eventT_B);
    maxtime=max(combotemp);
    maxtime=ceil(maxtime);%Times are in seconds
    timebins=0:0.001:maxtime;
    TbinsA=histc(spk_eventT,timebins);
    TbinsB=histc(spk_eventT_B,timebins);
    maxlag=200;
    [c1, plotlags]=xcorr(TbinsA,TbinsB,maxlag,'coeff');
    
    
    figure; plot(plotlags,c1);
    xlabel('Time Lag in ms', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    ylabel('Normalized Event Cross-Correlation', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    
    box off;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','inches');
    set(gcf,'Position',[1 1 10 7.5]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 7.5]);
    set(gcf, 'PaperPosition', [0 0 10 7.5]);
    set(gca,'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
    title([savename ' crosscorrelogram'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
    norm_histogram = gca;
    
end
%% plots spikes
figure;
if isempty(channel_B) ~= 1;
    subplot(6,1,[1 2]);
else
    subplot (3,1,[1 2]);
end

% let's remove super super huge artifacts for plotting: assume CC
idx = find(longV > .1);
longV(idx) = NaN;
idx = find(longV < -.2);
longV(idx) = NaN;

hold on
plot(longT, longV, 'k')
b = ylim;

yvect = [(b(1) + (.05*(b(2)-b(1))))  (b(1) + (.01*(b(2)-b(1))) )];
for i = 1:length(spk_eventT)
    plot([spk_eventT(i), spk_eventT(i)], yvect, 'r')
end


% plots drugline
if drugs.drugpresent == 1
    drugvect = [drugs.drugON drugs.drugOFF];
    yvaldrug = b(2) - (txtoffset_A*(b(2)-b(1)));
    yvect = [yvaldrug yvaldrug];
    plot(drugvect, yvect, 'b', 'LineWidth', 3)
    text(mean(drugvect), b(2) - ((txtoffset_A-.03)*(b(2)-b(1))), drugs.drugname, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
end
% plots second drug line
if drugs.drug2present == 1
    drugvect = [drugs.drugON2 drugs.drugOFF2];
    yvaldrug = b(2) - (txtoffset_B*(b(2)-b(1)));
    yvect = [yvaldrug yvaldrug];
    plot(drugvect, yvect, 'b', 'LineWidth', 3)
    text(mean(drugvect),b(2) -  ((txtoffset_B-.03)*(b(2)-b(1))), drugs.drugname2, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
end


% plots third drug line
if drugs.drug3present == 1
    drugvect3 = [drugs.drugON3 drugs.drugOFF3];
    yvaldrug = b(2) - (txtoffset_C*(b(2)-b(1)));
    yvect = [yvaldrug yvaldrug];
    plot(drugvect3, yvect, 'b', 'LineWidth', 3)
    % add text drug name
    text(mean(drugvect3), b(2) - (txtoffset_C*(b(2)-b(1))), drugs.drugname3, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
end



box off;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 10 7.5]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [15 12.5]);
set(gcf, 'PaperPosition', [0 0 10 7.5]);
if yOFF == 1
    set(gca, 'YTick', []);
    set(gca,'YColor','w');
    axis off;
end
% set(gca, 'XTick', 0:100:(longT(length(longV))+1), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
%xlabel('Time (s)');

xlim([0 longT(end)])

% plot pA scale bar

defxlim = xlim;
defylim = ylim;
scaleXx = [ (defxlim(2) - scaleX1) (defxlim(2))];
scaleXy = [ defylim(1) defylim(1)];
scaleYx = [ defxlim(2) defxlim(2)];
scaleYy = [ (defylim(1) + scaleY) (defylim(1))];

plot(scaleXx, scaleXy, 'k', 'Linewidth', 3)
plot(scaleYx, scaleYy, 'k', 'Linewidth', 3)

% yscale = (10^-3) * ymV_scalebar;
% yvect = [(mean(longV) - (yDOWN-3)*std(longV)) (mean(longV) - (yDOWN-3)*std(longV) + yscale)];
% plot([1, 1],yvect, 'k', 'linewidth',3)
% text(5, mean(yvect), [num2str(ymV_scalebar) ' mV' ], 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')

title([savename], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
hold off

%plots firing rate
if isempty(channel_B) ~= 1;
    subplot(6,1,3);
else
    subplot (3,1,3);
end
plot(binIwin, spk_bin_r, 'k')
set(gca, 'XTick', 0:round(longT(end)/10):longT(end), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
xlabel('Time (s)');
ylabel('Instantaneous Firing Rate (Hz)','Fontsize',10);
xlim([0 longT(end)])


% plots trace for B
if isempty(channel_B) ~= 1;
    subplot(6,1,[4 5]);
    
    plot(longT_B, longV_B, 'k'); hold on
    
    % plots drugline
    if drugs.drugpresent == 1
        drugvect = [drugs.drugON drugs.drugOFF];
        yvaldrug = mean(longV_B) + txtoffset*std(longV_B);
        yvect = [yvaldrug yvaldrug];
        plot(drugvect, yvect, 'b', 'LineWidth', 3)
        text(mean(drugvect), (yvaldrug + std(longV)), drugs.drugname, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
    end
    % plots second drug line
    if drugs.drug2present == 1
        drugvect = [drugs.drugON2 drugs.drugOFF2];
        yvaldrug = mean(longV_B) + txtoffset_B*std(longV_B);
        yvect = [yvaldrug yvaldrug];
        plot(drugvect, yvect, 'b', 'LineWidth', 3)
        text(mean(drugvect), (yvaldrug + std(longV)), drugs.drugname2, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
    end
    
    
    yvect = [(mean(longV_B) - 11*std(longV_B)) (mean(longV_B) - 9*std(longV_B))];
    for i = 1:length(spk_eventT_B)
        plot([spk_eventT_B(i), spk_eventT_B(i)], yvect, 'r')
    end
    axis off;
    ylim([(mean(longV_B) -13*std(longV_B)) (mean(longV_B) + 30*std(longV_B))])
    xlim([0 longT_B(end)])
    
    % plots firing rate for B
    subplot(6,1,6);
    plot(binIwinB, spk_bin_r_B, 'k')
    set(gca, 'XTick', 0:round(longT_B(end)/10):longT_B(end), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Instantaneous Firing Rate (Hz)','Fontsize',10);
    xlim([0 longT(end)])
    
    
    
    
end






IPSCraster = gca;
% dataPSC.pre_drugR = length(find(PSC_eventT < drugs.drugON))/(drugs.drugON); % PSC per second
% dataPSC.on_drugR = length(find(PSC_eventT > drugs.drugON & PSC_eventT < drugs.drugOFF))/(drugs.drugOFF - drugs.drugON); % PSC per second
%
if yes_save == 1
    saveas(IPSCraster, [savename ' Instantaneous firing rate.jpg'], 'jpg')
    %     save([savename, ' Inst_FR.mat'], 'longV', 'longT', 'spk_bin_ct', 'spk_bin_r', 'spk_eventT', 'reord_tracelist_V', 'Experiment', 'threshold', 'binT')
    save([savename, ' Inst_FR.mat'], '-regexp',  '^(?!Trace_.*$).')
    
end
savename
total_time = toc