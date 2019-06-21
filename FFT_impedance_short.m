% assume that there is only 1 sweep for each trace
% calculates impedance, phase angle, resonance frequency, resonance
% strength, and input resistance...
close all
clearvars -except Trace_*

% [FileName,PathName,FilterIndex] = uigetfile('*.mat');
%
% cd(PathName)
% load(FileName)

yes_save = 1;
graphtextoffset = 0;

TITLE = '190619 slice1 (1) AOB GC - 80mV exp .2-20Hz';

% which experiment
Experimentnum =     1    ;
Experiment = ['Trace_' num2str(Experimentnum) '_'];
% load wh:ich traces?
Trace = [  9:13    ];

% frequency cutoffs in ZAP
lowercutoff= .01; % .4
uppercutoff= 12; % 15
% lowercutoff= .2; % .4
% uppercutoff= 60; % 15

% % ZAP stimulus time in seconds and when it started in the trace
% ZAPtime =15;
% ZAPstart = 1; % time during trace when ZAP protocol started;
% currstepON = 17; % did current injection 12s into trace
% currsteptime = 1; % current injection lasted for 1s

ZAPstart = 5; % time during trace when ZAP protocol started;
ZAPtime = 20; % length of ZAP protocol;    10
currstepON = 1; % did current injection 12s into trace
currsteptime = 0.5; % current injection lasted for 1s

% ZAPstart = 5; % time during trace when ZAP protocol started;
% ZAPtime = 30; % length of ZAP protocol;    10
% currstepON = 1; % did current injection 12s into trace
% currsteptime = 1; % current injection lasted for 1s

target = .4 ;% Hz
% target = 2; % Hz
% target = lowercutoff;
% window size for filtering
dofilter = 1;
windowSize = 7;
resampF = 1;

% average and THEN do the FFT? 1 yes. 0 will do FFT and then average.
% average_THEN_FFT = 1;

% channel = 1 for voltage trace;   channel = 2 for current trace
% make list of voltage traces
% for i = 1:length(Trace)
%     reord_tracelist_voltage{i} = [Experiment, num2str(Trace(i)), '_', '1', '_', num2str(1)];
% end
% make list of current traces
for i = 1:length(Trace)
    reord_tracelist_current{i} = [Experiment, num2str(Trace(i)), '_', '1', '_', num2str(2)];
end

% calculate sampling rate
tmptmp = eval(reord_tracelist_current{1});
Fs = floor(length(tmptmp(:,1))/tmptmp(end,1))


%% Average traces and THEN perform FFT
% if average_THEN_FFT == 1;
for i = 1:length(reord_tracelist_current)
%     voltage = eval(reord_tracelist_voltage{i});
    current = eval(reord_tracelist_current{i});
%     toavg_volt(:,i) = voltage(:,2);
    toavg_curr(:,i) = current(:,2);
end

    prevlen = 0;

    for i = 1:length(Trace)
        tracelist = who([Experiment, num2str(Trace(i)), '_*_', num2str(1)]);
       for j = 1:length(tracelist)
                reord_tracelist{j} = [Experiment,num2str(Trace(i)), '_', num2str(j), '_', num2str(1)];
       end
        for j = 1:length(reord_tracelist)
            trace = eval(reord_tracelist{j});
            data_T = trace(:,1);
            
%%             % filter
            maxT = round(trace(end,1));
            Fs = length(trace(:,1))/maxT;
%             [bLP,aLP] = butter(2, 300/(Fs/2), 'low'); %300Hz LP filter 2nd order
%             trace(:,2) = (filtfilt(bLP,aLP, double(trace(:,2))));

%             data_V(:,i) = trace(:,2);
            
            toavg_volt(:, ( ((i-1) * prevlen) +j)) = trace(:,2);

%             plot(trace(:,1), trace(:,2));
            
            
        end
        prevlen = length(tracelist);
    end
    
avg_volt = mean(toavg_volt, 2);
avg_volt_mV = avg_volt.*10^3; % convert from V to mV
avg_curr = mean(toavg_curr, 2);
% filter
%     avg_volt_mV = filter(ones(1,50)/50,1,avg_volt_mV);

% crop out just the ZAP
toFFT_volt = avg_volt(((ZAPstart*Fs)+1):((ZAPstart+ZAPtime)*Fs));

toFFT_curr = avg_curr(((ZAPstart*Fs)+1):((ZAPstart+ZAPtime)*Fs));

% FFT to get impedance
numsamps=ZAPtime*Fs;
NFFT=2^nextpow2(numsamps);
freqs=(Fs/2)*linspace(0,1,NFFT/2+1);

[temp, highcutoff]=find(freqs>uppercutoff,1);
[temp, lowcutoff]=find(freqs>lowercutoff,1);
freqsnew=freqs(lowcutoff:highcutoff);
specsize=length(freqsnew);

%     ffV=fft(toFFT_volt, NFFT)/numsamps;
%     ffI=fft(toFFT_curr, NFFT)/numsamps;
ffV=fft(toFFT_volt)/numsamps;
ffI=fft(toFFT_curr)/numsamps;
%     currenth=hann(length(toFFT_curr));
%     voltageh=hann(length(toFFT_volt));

% currenth=nuttallwin(length(toFFT_curr));
% voltageh=nuttallwin(length(toFFT_volt));

%     ffV=fft(toFFT_volt.*voltageh, NFFT)/numsamps;
%     ffI=fft(toFFT_curr.*currenth, NFFT)/numsamps;

ffV2=ffV(lowcutoff:highcutoff);
ffI2=ffI(lowcutoff:highcutoff);

z1=ffV2./ffI2;
%     zmag=abs(z1)/ 10^6;
zmag = sqrt(imag(z1).^2 + real(z1).^2)/ 10^6;

phaseangle=atan(imag(z1)./real(z1));

% filtering/smoothing the numbers
if dofilter == 1;
    % Design the filter
    %         [bLP,aLP] = butter(6,100/(Fs/2), 'low'); %300Hz LP filter 2nd order
    %         impedance = (filtfilt(bLP,aLP,zmag)); %Low pass @300Hz
    %         phase = (filtfilt(bLP,aLP,phaseangle)); %Low pass @300Hz
    impedance = filter(ones(1,windowSize)/windowSize,1,zmag);
    phase = filter(ones(1,windowSize)/windowSize,1,phaseangle);
    
    
    %         impedance = filter(ones(1,windowSize)/windowSize,1,impedance);
    %         phase = filter(ones(1,windowSize)/windowSize,1,phase);
    %
    %         impedance = filter(ones(1,windowSize)/windowSize,1,impedance);
    %         phase = filter(ones(1,windowSize)/windowSize,1,phase);
    %         impedance = filter(ones(1,windowSize)/windowSize,1,impedance);
    %         phase = filter(ones(1,windowSize)/windowSize,1,phase);
    %
    impedance = resample(impedance, 1,resampF);
    phase = resample(phase, 1,resampF);
    freqsnew = resample(freqsnew, 1,resampF);
    
elseif dofilter == 0;
    impedance = zmag;
    phase = phaseangle;
end

% get resonance freq for cell
%     newidx_imp = find(freqsnew > 1 & freqsnew < 7);
%         newidx_imp = find(freqsnew > 1);
newidx_imp = find(freqsnew > target);
maximped_val = max(impedance(newidx_imp));
%     maximped_val = max(impedance);
idx = find(impedance == maximped_val);
maximped_freq = freqsnew(idx);

% calculate Q (ratio of impedance @ resonance freq over impedance @ lowest
% frequency (.5Hz))
%     target = .01; % Hz
idx = find((target-.05) < freqsnew & freqsnew < (target +.05));

% temp = abs(freqsnew(idx)-target);
% [idx idx] = min(temp);

min_idx = min(impedance(idx));
idx = find(impedance == min_idx);
lowfreqimp_val = impedance(idx);
Qfactor = maximped_val/lowfreqimp_val;

freqsnew_forplotting = freqsnew(idx:length(freqsnew));
impedance_forplotting = impedance(idx:length(impedance));
phase_forplotting = phase(idx:length(phase));




% get time(s) for xaxis for plotting
temp = eval(reord_tracelist_current{1});
xtime = temp(:,1);

% calculate input resistance --> average .25s of data from before current
% injection on and before current injection off
avg_baseline_volt = mode([avg_volt(((Fs*(currstepON - .25))-1):(((Fs*currstepON)-1)))]);
avg_curr_inj_volt = mode([avg_volt(((Fs*((currstepON + currsteptime)-.25)-1)):((Fs*(currstepON + currsteptime)-1)))]);

avg_baseline_curr = mode([avg_curr(((Fs*(currstepON - .25))-1):(((Fs*currstepON-1))))]);
avg_curr_inj_curr = mode([avg_curr(((Fs*((currstepON + currsteptime)-.25)-1)):((Fs*(currstepON + currsteptime)-1)))]);

deltaV = avg_curr_inj_volt - avg_baseline_volt;
deltaC = avg_curr_inj_curr - avg_baseline_curr;
inputR = abs((deltaV/deltaC))/10^9; % input resistance in giga Ohms



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalize impedance profile




%% make a plot

IMPOTENCE = figure;

% subplot(3,1,1); plot(xtime, avg_volt_mV, 'Color', 'k', 'LineWidth', 1.5)
subplot(2,1,1); plot(xtime, avg_volt_mV, 'Color', 'k', 'LineWidth', 1.5)
ylabel('Membrane Potential (mV)',  'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time',  'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
yup = abs(max(avg_volt_mV)/10) + max(avg_volt_mV);
ydown = min(avg_volt_mV) - abs(max(avg_volt_mV)/10);
ylim([ydown yup])
text(xtime(((currstepON)*Fs)), min(avg_volt_mV), ['Input Resistance = ', num2str(inputR), 'G\Omega'], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment', 'left')
title('Impedance Profile', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

freqsnew_forplotting = freqsnew_forplotting';


% subplot(3,1,2); plot(freqsnew_forplotting,impedance_forplotting, 'Color', 'k', 'LineWidth', 1.5);
subplot(2,1,2); plot(freqsnew_forplotting,impedance_forplotting, 'Color', 'k', 'LineWidth', 1.5);
ylabel('Impedence (M\Omega)',  'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Frequency (Hz)',  'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
xlim([lowercutoff uppercutoff])
text((freqsnew(floor(length(freqsnew)*.4))), (maximped_val+graphtextoffset), ['Resonance Freq = ', num2str(maximped_freq), 'Hz'], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment','top', 'HorizontalAlignment', 'center');
text((freqsnew(floor(length(freqsnew)*.8))), (maximped_val+graphtextoffset), ['Q factor = ', num2str(Qfactor)], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'VerticalAlignment','top', 'HorizontalAlignment', 'center');

% subplot(3,1,3); plot(freqsnew_forplotting,phase_forplotting , 'Color', 'k', 'LineWidth', 1.5);
% ylabel('Phase angle (radians)', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('Frequency (Hz)',  'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
% xlim([lowercutoff uppercutoff])


% IMPOTENCE = figure(1);
if yes_save == 1
    saveas(IMPOTENCE, [TITLE '.jpg'], 'jpg');
    save([TITLE, '.mat'], 'impedance', 'target', 'phase','freqsnew','inputR', 'maximped_freq', 'maximped_val', 'toavg_volt', 'reord_tracelist_current','lowfreqimp_val', 'Qfactor', 'avg_volt', 'avg_volt_mV', 'avg_curr', 'xtime', 'windowSize', 'resampF',  'freqsnew_forplotting', 'impedance_forplotting', 'phase_forplotting');

%     save([TITLE, '.mat'], 'impedance', 'target', 'phase','freqsnew','inputR', 'maximped_freq', 'maximped_val', 'toavg_volt', 'reord_tracelist_current','lowfreqimp_val', 'Qfactor', 'avg_volt', 'avg_volt_mV', 'avg_curr', 'xtime', 'windowSize', 'resampF', 'to_avg', 'freqsnew_forplotting', 'impedance_forplotting', 'phase_forplotting');
end