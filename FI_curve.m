close all
clearvars -except Trace_*

tic
%% settings and shit
savename = '190619 slice1 (1) AOB GC';
yes_save = 1;

% which experiment
Experimentnum =         1          ;   % 
% load which trace?
Trace =                 5          ;   % 
channel =               1          ;
currchannel =           2          ;

yOFF = 1;

% Threshold for spikes
threshold = -.0000;
spk.slopethresh = .01; % .1; 3 for Fs?
inst_spkr_bin = .1; % 10 ms instantaneous spike rate bins

% when the current step starts, when it` ends
stimstart = 1; % seconds
stimend = 3;

%do sag? if it has Ih, calculate the sag amplitude
dosag = 1;
sagonly = 0;

% plot taus?
plottau = 1;
taustop = .25;

%%  set up the list of traces and find Fs
Experiment = ['Trace_' num2str(Experimentnum) '_'];
tracelist = who([Experiment, num2str(Trace), '_*_', num2str(channel)]);

for i = 1:length(tracelist)
    reord_tracelist_V{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(channel)];
    reord_tracelist_I{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(currchannel)];
end

% find Fs
tmp = eval(reord_tracelist_V{1});
maxT = tmp(end,1); % time index for end
Fs = floor(length(tmp(:,1))/maxT);

%% detect spikes and get current values!!!
for i = 1:length(reord_tracelist_V)
    % spike detection
    tmpV = eval(reord_tracelist_V{i});
    
    data_V{i} = tmpV(:,2);
    data_T{i} = tmpV(:,1);
    
    [pks, lct] = findpeaks(tmpV(:,2), 'MinPeakHeight', threshold, 'MinPeakDistance', round(.001*Fs));
    spk_peventsT = tmpV(lct,1);
    spk.T{i,:} = num2cell(spk_peventsT);
    spk.height{i,:} = num2cell(pks);
    
    % get current injection
    tmpI = eval(reord_tracelist_I{i});
    baseI = mean(tmpI(1:(Fs*stimstart), 2));
    stimI = mean(tmpI( Fs*stimstart: Fs*stimend, 2));
    spk.stepI(i) = (stimI - baseI) * 10^12; % converts to pA
    
    % get cellprop.sag amplitude
    if spk.stepI(i) < -1
        cellprop.Vmin(i) = min(tmpV( (Fs*stimstart): (Fs*stimend), 2));
        cellprop.Vsteady(i) = mean(tmpV( (Fs* (stimend - .25)) : (Fs*stimend), 2));
        
        cellprop.Vbase(i) = mean(tmpV( 1 : Fs*stimstart, 2));
        
        cellprop.sag(i) = cellprop.Vmin(i) - cellprop.Vsteady(i);
        cellprop.sag(i) = cellprop.sag(i) * 1000; % converts to mV
        
        % calculates input resistance
        cellprop.inputR(i) = ( (cellprop.Vsteady(i) - cellprop.Vbase(i) ) / spk.stepI(i)) * 10^3; % IN GIGA OHMS
        
        % calculate tau
        xtofit = tmpV( (Fs*stimstart): (Fs*(stimstart+taustop)), 1);
        ytofit = tmpV( (Fs*stimstart): (Fs*(stimstart+taustop)), 2);
        xmod = xtofit - xtofit(1);
        ymod = ytofit - ytofit(end);
        
        % double exp
        %         ft = fittype( 'a*exp(-x/b)+c*exp(-x/d)', 'independent', 'x', 'dependent', 'y' );
        %         taufit = fit(xmod,ymod, ft);
        %         tauw(i) = ( taufit.b * (taufit.a/(taufit.a+taufit.c)) ) + ( taufit.d * (taufit.c/(taufit.a+taufit.c)));
        % %         tauw(i) = 1/tauw(i);
        %         tmp = tauw(i);
        %         tmp
        
        % single exp
        ft = fittype( 'a*exp(-x*b)', 'independent', 'x', 'dependent', 'y' );
        taufit = fit(xmod,ymod, ft);
        %         tauw(i) = ( taufit.b * (taufit.a/(taufit.a+taufit.c)) ) + ( taufit.d * (taufit.c/(taufit.a+taufit.c)));
        %         tauw(i) = 1/tauw(i);
        tau(i) = taufit.b;
        
        if plottau == 1
            figure;plot(taufit, xmod,ymod); hold on; title(['tau = ' num2str(tau(i)) ' . cellprop.Vmin = ' num2str(cellprop.Vmin(i)) '  .  cellprop.Vbase = ' num2str(cellprop.Vbase(i))])
            
            prompt = input('Keep this fit? (1 = yes, 0 = no): ');
            if prompt == 0
                tau(i) = NaN;
                
                
            end
            if ~isnan(tau(i))
                cellprop.Capac(i) = tau(i)/cellprop.inputR(i)  * 10^-6;
            end
            close all
        end
        
    end
    
end

if sagonly == 0
%% make PSTH
% binT = .1; % .1 s bin size
binT = inst_spkr_bin;
spk.binIwin = [binT:binT:ceil(maxT/binT)*binT]; % or you could floor it?

spk.bins = zeros(length(reord_tracelist_V), length(spk.binIwin));

for i = 1:length(reord_tracelist_V)
    tmp = cell2mat(spk.T{i,:});
    if isempty(tmp) ~=1
        spk.bins(i,:) = histc(tmp, spk.binIwin);
        spk.binsrate(i,:) = spk.bins(i,:)/binT;
    end
end

PSTH = sum(spk.bins, 1);
%% Make F-I plot data!!!
for i = 1:length(reord_tracelist_V)
    tmp = cell2mat(spk.T{i,:});
    idx = find(tmp > stimstart & tmp < stimend);
    
    if isempty(idx) == 1
        spkwincount(i) = 0;
        spk.spkwinrate(i) = 0;
        spk.CVisi(i) = 0;
    else
        spkstimwin = tmp(idx);
        spkwincount(i) = length(spkstimwin);
        spk.spkwinrate(i) = spkwincount(i)/(stimend-stimstart);
        spk.ISI{i} = diff(spkstimwin);
        spk.CVisi(i) = std(spk.ISI{i})/mean(spk.ISI{i});
        %         spk.CVisi(i) = ( max(spk.ISI{i}) - min(spk.ISI{i}) ) /mean(spk.ISI{i});
        %         spk.CVisi(i) = median(spk.ISI{i}) /mean(spk.ISI{i});
    end
    
end


%% calculate spk.slopes bro!!!!
% Method A
% spk.slopethresh = .25;
% dspk.stepI = diff(spk.stepI);
% dspkr = diff(spk.spkwinrate);
%
% start = dspkr(1)/dspk.stepI(1);
% spk.slopes = [];
% spk.slopeidx = [];
% for i = 2:length(dspkr)
%     tmp = dspkr(i)/dspk.stepI(i);
%     if tmp < (1+spk.slopethresh)*start && tmp > (1-spk.slopethresh)*start
%         spk.slopes = [spk.slopes tmp];
%         spk.slopeidx = [spk.slopeidx i];
%     end
%     start = tmp;
% end
% spk.slopes
% spk.slope_avg = mean(spk.slopes) % hz
%% Rich's and Larry Abbott's method from Chance et al 2002
dspkr = diff(spk.spkwinrate);
idx = find(dspkr > spk.slopethresh);
newidx = [idx (idx(end)+1)];
xvectfit = spk.stepI(newidx);
yvectfit = spk.spkwinrate(newidx);
P = polyfit(xvectfit, yvectfit,1);
yfit = P(1)*xvectfit+P(2); % equation for straight line
spk.slope = P(1);
spk.intercept = P(2);

%% FIGURE TIME

figure;
subplot(5,2,[ 1:4]);
hold on
title([savename ' FI'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');

pcolorsss = colormap(jet(length(reord_tracelist_V)));
for i = 1:length(reord_tracelist_V)
    tmp = eval(reord_tracelist_V{i});
    plot(tmp(:,1), tmp(:,2), 'Color', pcolorsss(i,:))
end
xlabel('Time (s)','FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
ylabel('V','FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
xlim([0 maxT])
set(gca, 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');

subplot(5,2,[5 6]);
yvect = [0 .8];
hold on
for i = 1:length(reord_tracelist_V)
    tmp = cell2mat(spk.T{i,:});
    for j = 1:length(tmp)
        plot([tmp(j) tmp(j)], (yvect - i), 'Color', pcolorsss(i,:)) % subtracting makes it so that first trace is top
    end
end

if yOFF == 1
    set(gca, 'YTick', []);
    set(gca,'YColor','w');
end
set(gca, 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
xlabel('Time (s)');
xlim([0 maxT])

% plots PSTH
subplot (5,2,[7 8]);
hold on
for i = 1:length(spk.binsrate(:,1))
    tmp = cell2mat(spk.T{i,:});
    if isempty(tmp) ~=1
        plot(spk.binIwin, spk.binsrate(i,:), 'Color', pcolorsss(i,:))
    end
end
% plot(spk.binIwin, PSTH, 'k')
set(gca, 'XTick', 0:1:maxT, 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
xlabel('Time (s)');
ylabel('Inst Spike Rate (Hz)','FontSize',10);
xlim([0 maxT])
text(  mean(spk.binIwin), mean(ylim), ['Spike rate bin size = ' num2str(binT) 's'], 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

% RASTER = gca;
% hold off

% %% plot FI curve!
% figure
subplot(5,2,[9]);
hold on

% title([savename ' FI curve'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');

plot(spk.stepI, spk.spkwinrate, '-ok')
xlabel('Current injection (pA)', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
ylabel('Spike rate (Hz)', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
text(min(spk.stepI), mean(spk.spkwinrate), ['spk.slope = ' num2str(spk.slope) ' Hz/pA'], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment', 'left')
plot(xvectfit, yfit, 'b-.', 'Linewidth', 2)
plot(xvectfit, yvectfit, 'bo', 'Markersize', 10)

% 
% subplot(3,1,2);
% plot(spk.stepI, spkwincount, '-ok')
% xlabel('Time (s)');
% ylabel('Number of spikes', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
% xlabel('Current injection (pA)', 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
% set(gca, 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');

subplot(5,2,[10]);
plot(spk.stepI, spk.CVisi, '-ok')
xlabel('Time (s)');
ylabel('CV_I_S_I', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
xlabel('Current injection (pA)', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');


box off;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 11 9.5]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [11 9.5]);
set(gcf, 'PaperPosition', [0 0 11 7.5]);
set(gca, 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');

FIcurve = gca;

end
%% 
if sagonly ==1
    
figure;
% subplot(5,2,[ 1:4]);
hold on
title([savename ' FI'], 'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');

pcolorsss = colormap(jet(length(reord_tracelist_V)));
for i = 1:length(reord_tracelist_V)
    tmp = eval(reord_tracelist_V{i});
    plot(tmp(:,1), tmp(:,2), 'Color', pcolorsss(i,:))
end
xlabel('Time (s)','FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
ylabel('V','FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
xlim([0 maxT])
set(gca, 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
FIcurve = gca;

end
%%


if yes_save == 1
%     saveas(RASTER, [savename ' CCIV PSTH.jpg'], 'jpg')
    saveas(FIcurve, [savename ' CCIV FI curve.jpg'], 'jpg')
    save([savename, ' CCIV.mat'], '-regexp',  '^(?!Trace_.*$).')
end
% cellprop.sag
total_time = toc;