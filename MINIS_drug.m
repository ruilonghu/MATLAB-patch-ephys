close all
clearvars -except Trace_*
warning('off')
tic

%% settings and shit
savename =  '190619 slice1 (1) AOB GC - EPSC' ;
yes_save =                    1               ;

% which experiment
Experimentnum =               1               ;
% load which trace?
Trace =  [                    4              ]; %15 16 17
channel =                     1               ;

% are you doing EPSCs?
EPSCs =                       1               ; % 0 for IPSCs, 1 for EPSCs, 2 for KCl
EPSP =                        0               ;

if EPSP == 0
    if EPSCs == 1
        base_threshold = (10^-12) *     3         ; % 2to10 pA for EPSCs
        high_threshold = (10^-12) *     2.5         ; % make it 1-2pA smaller than base_threshold
    elseif EPSCs == 0
        base_threshold = (10^-12) *     8         ; % 5to10 pA for CsG IPSCs
        high_threshold = (10^-12) *     6         ; % make it 1-2pA smaller than base_threshold
    elseif EPSCs == 2
        base_threshold = (10^-12) *     12        ; % 8to20 pA for CsClQX IPSCs
        high_threshold = (10^-12) *     10       ; % make it 1-2pA smaller than base_threshold
    end
elseif EPSP == 1
    base_threshold = (10^-3) * 1; % for EPSPs
    high_threshold = (10^-3) * .1 ; % 10pA threshold to count IPSCs in the superfiltered trace - 2to5 pA for EPSCs
end

%% drug info
drugs.drugpresent = 0;
drugs.drugname = 'muscarine 10uM';
drugs.drugON = 60 * 10; % seconds (from minutes)
drugs.drugOFF = 60 * 15; % seconds
%
% %
% drugs.drugpresent = 0;
% drugs.drugname = 'clonidine 10 uM';
% drugs.drugON = 60 * 2; % seconds (from minutes)
% drugs.drugOFF = 60 * 10; % seconds
%
% drugs.drugpresent = 0;
% drugs.drugname = 'TTX  0.5 uM';
% drugs.drugON =  60 * 3; % seconds (from minutes)
% drugs.drugOFF = 60 * 7; % seconds

% calculate drug xxx seconds past drug turns on
drugs.drugonDELAY = 60 * 5 ; % how long to wait before start analyzing for drug effect
drugs.drugonANAL = 60 * 3; % amount of time window used to analyze drug effect

drugs.analwin1_pre = [(drugs.drugON - drugs.drugonANAL) , drugs.drugON ];
drugs.analwin1_post = [(drugs.drugON + drugs.drugonDELAY) , (drugs.drugON + drugs.drugonDELAY + drugs.drugonANAL)];

drugs.drug2present = 0;
drugs.drugname2 = '4DAMP 100nM';
drugs.drugON2 = 60 * 25; % seconds (from minutes)
drugs.drugOFF2 = 60 * 43; % seconds

drugs.drug3present = 0;
drugs.drugname3 = 'muscarine 10uM';
drugs.drugON3 = 60 * 36; % seconds (from minutes)
drugs.drugOFF3 = 60 * 40; % seconds

% drug text offset as a fraction of the full y axis
txtoffset_A = .1; % default .1
txtoffset_B = .15;
txtoffset_C = .2;
%% other params
% do autocorrelation?
autoCORR = 0;
msbin = 0; % using seconds or milliseconds

% do wave properties?
do_wave = 1;

% yaxis off?
yOFF =      0    ;
% scale bar
y_scalebar = 20; % pA for scale bar
yes_scalebar = 1;

%%  set up the list of traces and find Fs
Experiment = ['Trace_' num2str(Experimentnum) '_'];
% tracelist = who([Experiment, num2str(Trace), '_*_', num2str(channel)]);
% %
% for i = 1:length(tracelist)
%     reord_tracelist_I{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(channel)];
% end

% tmp = eval(reord_tracelist_I{1});
% maxT = tmp(end,1); % time index for end
% Fs = length(tmp(:,1))/maxT;

tmp = eval([Experiment, num2str(Trace(1)), '_1_', num2str(1)]);
maxT = round(tmp(end,1));
Fs = length(tmp(:,1))/maxT;

%% stack into long trace and do stuff

data_T = [];
data_I = [];
tmpdata_T = [];
tmpdata_V = [];
lastTi = 0;
for i = 1:length(Trace)
    tracelist = who([Experiment, num2str(Trace(i)), '_*_', num2str(channel)]);
    
    longT = [];
    longV = [];
    
    lastTj = 0;
    for j = 1:length(tracelist)
        tmp = eval([Experiment, num2str(Trace(i)), '_', num2str(j), '_', num2str(channel)]);
        tmpdata_T = tmp(:,1);
        tmpdata_V = tmp(:,2);
        
        longT = vertcat(longT, (tmpdata_T + lastTj));
        longV = vertcat(longV, tmpdata_V) ;
        
        lastTj = longT(end); % it'll be the end time for hte next loop iteration
        
    end
    
    data_T = vertcat(data_T, (longT + lastTi));
    data_I = vertcat(data_I, longV);
    
    lastTi = data_T(end);
end

longT = data_T;
longI_p = data_I;
%% detect PSCs

[PSC_event, PSCidx, normI, longI, wavemat, dataPSC, eventTvect] = detectPSCs(longI_p, longT, base_threshold, high_threshold, Fs, EPSCs);

wavemat = wavemat';
%% do things like look at amplitudes, rise times, and decays

dataPSC.amp_avg = nanmean(dataPSC.amplitude);
dataPSC.amp_std = nanstd(dataPSC.amplitude);
dataPSC.riseT_avg = nanmean(dataPSC.riseT);
dataPSC.riseT_std = nanstd(dataPSC.riseT);
dataPSC.decayT_avg = nanmean(dataPSC.decayT);
dataPSC.decayT_std = nanstd(dataPSC.decayT);
dataPSC.tau_avg = nanmean(dataPSC.tau);
dataPSC.tau_std = nanstd(dataPSC.tau);


%% calculate some numbers if drug is present

if drugs.drugpresent == 1
    % calculate drug times too
    
    % get indices for pre drug calculation and on drug effect calculation
    drugs.preidx = find(PSC_event > drugs.analwin1_pre(1)  & PSC_event < drugs.analwin1_pre(2));
    drugs.postidx =  find(PSC_event > drugs.analwin1_post(1) & PSC_event < drugs.analwin1_post(2));
    %
    
    dataPSC.pre_drugFREQ = length(drugs.preidx)/drugs.drugonANAL; % PSC per second
    
    if do_wave == 1
        dataPSC.pre_amp = dataPSC.amplitude(drugs.preidx);
        dataPSC.pre_amp_avg = nanmean(dataPSC.pre_amp);
        dataPSC.pre_riseT = dataPSC.riseT(drugs.preidx);
        dataPSC.pre_riseT_avg = nanmean(dataPSC.pre_riseT);
        dataPSC.pre_decayT = dataPSC.decayT(drugs.preidx);
        dataPSC.pre_decayT_avg = nanmean(dataPSC.pre_decayT);
        dataPSC.pre_tau = dataPSC.tau(drugs.preidx);
        dataPSC.pre_tau_avg = nanmean(dataPSC.pre_tau);
    end
    
    
    dataPSC.on_drugFREQ = length(drugs.postidx)/drugs.drugonANAL; % PSC per second
    
    if do_wave == 1
        dataPSC.post_amp = dataPSC.amplitude(drugs.postidx);
        dataPSC.post_amp_avg = nanmean(dataPSC.post_amp);
        
        dataPSC.post_riseT = dataPSC.riseT(drugs.postidx);
        dataPSC.post_riseT_avg = nanmean(dataPSC.post_riseT);
        
        dataPSC.post_decayT = dataPSC.decayT(drugs.postidx);
        dataPSC.post_decayT_avg = nanmean(dataPSC.post_decayT);
        
        dataPSC.post_tau = dataPSC.tau(drugs.postidx);
        dataPSC.post_tau_avg = nanmean(dataPSC.post_tau);
        
        
        [h p] = ttest2(dataPSC.pre_amp, dataPSC.post_amp);
        dataPSC.significance.amp = [h p];
        [h p] = ttest2(dataPSC.pre_riseT, dataPSC.post_riseT);
        dataPSC.significance.riseT = [h p];
        [h p] = ttest2(dataPSC.pre_decayT, dataPSC.post_decayT);
        dataPSC.significance.decayT = [h p]
        %
        [h p] = ttest2(dataPSC.pre_tau, dataPSC.post_tau);
        dataPSC.significance.tau = [h p];
        
        %     if EPSP == 1
        %         [h p] = ttest2(dataPSC.pre_halfwidth, dataPSC.post_halfwidth);
        %         dataPSC.significance.halfwidth = [h p];
        %     end
        
        dataPSC.nancount = sum(isnan(dataPSC.riseT));
        dataPSC.nanratio = sum(isnan(dataPSC.riseT))/length(dataPSC.riseT);
        dataPSC
        dataPSC.significance
        savename
    end
else
    dataPSC
    savename
end


%% do autocorrelation
if autoCORR == 1
    
    if drugs.drugpresent == 0
        combotemp=PSC_event;
        maxtime=max(combotemp);
        maxtime=ceil(maxtime);%Times are in seconds
        
        if msbin == 1
            tbinwin = .01;
            maxlag = 5000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins= 0:tbinwin:maxtime;
        TbinsA=histc(PSC_event,timebins);
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
        combotemp = PSC_event(drugs.preidx);
        
        if msbin == 1
            tbinwin = .001;
            maxlag = 5000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins=  drugs.analwin1_pre(1):tbinwin:drugs.analwin1_pre(2);
        TbinsA=histc(combotemp,timebins);
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
        
        combotemp = PSC_event(drugs.postidx);
        
        if msbin == 1
            tbinwin = .001;
            maxlag = 5000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins=  drugs.analwin1_post(1):tbinwin:drugs.analwin1_post(2);
        TbinsA=histc(combotemp,timebins);
        
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
        
        drugs.postidx = find(PSC_event > (longT(end) - 180) & PSC_event < (longT(end)));
        
        combotemp = drugs.postidx;
        
        maxtime=max(combotemp);
        maxtime=ceil(maxtime);%Times are in seconds
        
        if msbin == 1
            tbinwin = .001;
            maxlag = 5000; %ms
        else
            tbinwin = 1;
            maxlag = 60; %s
        end
        
        timebins=0:tbinwin:maxtime;
        TbinsA=histc(PSC_event,timebins);
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
%% Instantaneous PSC event rate histogram

combotemp = PSC_event;
maxtime=max(combotemp);
maxtime=ceil(longT(end));%Times are in seconds

InstTbin = 1; % this allows for 1 s bins

timebins=0:InstTbin:maxtime;
PSCbins=histcounts(PSC_event,timebins)'; % this is basically instantaneous event rate.

%% plots IPSCs
% figure; plot(longT, normI, 'k')
figure;
subplot(3,1,1:2)
plot(longT, longI, 'k')

hold on
scatter(PSC_event, longI(PSCidx), 'o', 'r')

b = ylim;

% plots drug line
if drugs.drugpresent == 1
    drugvect = [drugs.drugON drugs.drugOFF];
    yvaldrug = b(2) - (txtoffset_A*(b(2)-b(1)));
    yvect = [yvaldrug yvaldrug];
    plot(drugvect, yvect, 'b', 'LineWidth', 3)
    % add text drug name
    text(mean(drugvect), b(2) - (txtoffset_A*(b(2)-b(1))), drugs.drugname, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
end

% plots second drug line
if drugs.drug2present == 1
    drugvect2 = [drugs.drugON2 drugs.drugOFF2];
    yvaldrug = b(2) - (txtoffset_B*(b(2)-b(1)));
    yvect = [yvaldrug yvaldrug];
    plot(drugvect2, yvect, 'b', 'LineWidth', 3)
    % add text drug name
    text(mean(drugvect2), b(2) - (txtoffset_B*(b(2)-b(1))), drugs.drugname2, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold')
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
set(gcf, 'PaperSize', [10 7.5]);
set(gcf, 'PaperPosition', [0 0 10 7.5]);
if yOFF == 1
    set(gca, 'YTick', []);
    set(gca,'YColor','w');
end

set(gca, 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');

xlabel('Time (s)');
% general axis limits
xlim([0 longT(end)])

if yes_scalebar == 1
    % plot pA scale bar
    if EPSP == 0
        yscale = (10^-12) * y_scalebar;
    elseif EPSP == 1
        yscale = (10^-3) *  y_scalebar;
    end
    
    yvect = [(min(longI) - 4*std(longI))     ( (min(longI) - 4*std(longI) + yscale))];
    
    plot([1, 1],yvect, 'k', 'linewidth',3)
    if EPSP == 0
        text(5, mean(yvect), [num2str(y_scalebar) ' pA' ], 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
    elseif EPSP == 1
        text(5, mean(yvect), [num2str(y_scalebar) ' mV' ], 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold')
    end
end

title([savename], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');

ttbins = timebins(2:end)';
% plot instantaneous event rate
subplot(3,1,3)
plot(ttbins, PSCbins, 'k')
set(gca, 'XTick', 0:round(longT(end)/15):(longT(length(normI))+1), 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
xlim([0 longT(end)])
xlabel('Time (s)');
ylabel('Hz');

IPSCraster = gca;

hold off

%%
if yes_save == 1
    saveas(IPSCraster, [savename ' raster.jpg'], 'jpg')
    save([savename '.mat'], 'wavemat', 'dataPSC', 'base_threshold', 'high_threshold', 'EPSCs', 'longT', 'longI', 'normI','PSC_event', 'tracelist', 'Experiment', 'dataPSC', 'drugs', 'Fs','eventTvect','ttbins', 'PSCbins')
    
end

total_time = toc