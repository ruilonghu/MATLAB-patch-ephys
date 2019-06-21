function [PSC_event, PSCidx, normI_sm, longI, wavemat, dataPSC, eventTvect] = detectPSCs(longI_p, longT, base_threshold, high_threshold, Fs, EPSCs);
%Detects PSCs Summary of this function goes here
%   Detailed explanation goes here

% threshold = -threshold;

%% pass long trace through 300Hz low bass Butter filter
% if doFILTER == 1;
%     [bLP,aLP] = butter(2, 200/(Fs/2),'low'); %300Hz LP filter 2nd order
%     eval(sprintf('longI = (filtfilt(bLP,aLP,longI_p));')); %Low pass @300Hz
% else
%     longI = longI_p;
% end
% keyboard
longI = longI_p;

if EPSCs == 0
    %     [bLP,aLP] = butter(2, 10/(Fs/2),'low'); %300Hz LP filter 2nd order
    %         [bLP,aLP] = butter(2, 60/(Fs/2), 'low'); %300Hz LP filter 2nd order
    %     [bLP,aLP] = butter(2, 600/(Fs/2), 'low'); %300Hz LP filter 2nd order  ;  500 for FS epscs
    [bLP,aLP] = butter(2, 500/(Fs/2), 'low'); %500Hz LP filter 2nd order  ;  500 for FS epscs
    
elseif EPSCs == 2
    [bLP,aLP] = butter(2, 500/(Fs/2),'low'); %500Hz LP filter 2nd order
elseif EPSCs == 1
    [bLP,aLP] = butter(2, 500/(Fs/2), 'low'); %500Hz LP filter 2nd order  ;  500 for FS epscs
    %         [bLP,aLP] = butter(2, 1000/(Fs/2), 'low'); %500Hz LP filter 2nd order  ;  500 for FS epscs
    
    %     [bLP,aLP] = butter(2, 2/(Fs/2),'low'); %300Hz LP filter 2nd order
    %     [bLP,aLP] = butter(2, 60/(Fs/2),'low'); %300Hz LP filter 2nd order
    %     [bLP,aLP] = butter(2, 20/(Fs/2),'low'); %300Hz LP filter 2nd order
    
end
% default butter(2, 200/(Fs/2), 'low');

%% super lowpass time re-adjust
slptreadjust = 0;

%%

longI_sm = filtfilt(bLP,aLP,longI);  %Low pass @100Hz FOR SUPER FILTER
% longI_sm = smooth(longI_sm,'rlowess');

if EPSCs == 0
    longI_inv_sm = longI_sm;
    longI_inv = longI;
else
    longI_inv_sm = -longI_sm;% inverts so peaks are positive
    longI_inv = -longI;
end


%% find baseline values in 10s bins and then find IPSCs
binT = Fs * 1; % .5 s bin size
binIwin = [1:binT:length(longT)];
for i = 1:length(binIwin)
    %find baseline
    if i == length(binIwin)
        idxstop = length(longT);
    else
        idxstop = binIwin(i+1);
    end
    tmpI = longI_inv_sm(binIwin(i):idxstop);
    tmpT = longT(binIwin(i):idxstop);
    maxI = max(tmpI);
    minI = min(tmpI);
    histrange = [minI:(maxI - minI)/1000:maxI];
    Ihist = histc(tmpI, histrange);
    [val idx] = max(Ihist);
    if isempty(idx)
        keyboard
    end
    baseI = histrange(idx(1));
    % normalize baseline
    normI_sm(binIwin(i):idxstop) = tmpI - baseI;
    
    %find baseline and redo for the non-smoothed
    if i == length(binIwin)
        idxstop = length(longT);
    else
        idxstop = binIwin(i+1);
    end
    tmpI = longI(binIwin(i):idxstop);
    tmpT = longT(binIwin(i):idxstop);
    maxI = max(tmpI);
    minI = min(tmpI);
    histrange = [minI:(maxI - minI)/1000:maxI];
    Ihist= histc(tmpI, histrange);
    [val idx] = max(Ihist);
    baseI = histrange(idx(1));
    
    % normalize baseline
    normI(binIwin(i):idxstop) = tmpI - baseI;
end


% keyboard

if EPSCs ~= 0
    normI = -normI;
end

%% turn those thresholded values into single PSC events

% if EPSCs == 1
%     mindist = Fs* .0035;
% elseif EPSCs == 0 % for outward IPSCs
%     mindist = Fs* .005; % default .01?
% elseif EPSCs == 2 % for KCl inwards
%     mindist = Fs* .005;
% end


% [pks, PSCidx] = findpeaks(normI_sm, 'MinPeakHeight', threshold, 'MinPeakDistance', mindist);
% [pks, PSCidx] = findpeaks(normI_sm, 'MinPeakHeight', base_threshold, 'MinPeakProminence', high_threshold, 'MinPeakDistance', mindist);
[pks, PSCidx] = findpeaks(normI_sm, 'MinPeakHeight', base_threshold, 'MinPeakProminence', high_threshold);

PSC_event_tmp = longT(PSCidx);
% PSC_event_new = PSC_event;
%%
%% let's try the matrix deconv way
%

%% use average to generate filter
% % [pks, PSCidx] = findpeaks(longI_inv(1:), 'MinPeakHeight', threshold, 'MinPeakDistance', mindist);
% [pks, PSCidx] = findpeaks(longI_inv((Fs*2*60):(Fs*3*60)), 'MinPeakHeight', threshold, 'MinPeakDistance', mindist);
%
% PSC_event = longT(PSCidx);
%
% window1=[0.01 0.05]; %First value is for the part before, and the second value is for the values after the threshold crossing
% numevents=length(PSC_event);
% window2=round(window1*Fs);
% windowlen=sum(window2) + 1;
% %
% tmpwavemat=zeros(numevents,windowlen);
% % figure; hold on
% for ik=1:numevents
%     i1=PSC_event(ik)*Fs;
%     b1=i1-window2(1);
%     b2=i1+window2(2);
%     if b1>0 && b2<length(longI_inv)
%         tempwav=longI_inv(b1:(b2));
%         tmpwavemat(ik,:) = tempwav;
%         %    plot(wavemat(ik,:))
%     end
%
% end
% %average and normalize the templatewave
% templatewave = mean(tmpwavemat,1);
% templatewave = -templatewave/max(templatewave);
% templatewave = templatewave - max(templatewave);
% templatewave = -templatewave/min(templatewave);

% keyboard
% generate filter using means
% avgwave = fliplr(mean(wavemat, 1));

% plot(avgwave)

%% make your own damn filter
% f_len = round(Fs * .03);
% f_rise = round(Fs * .002); % 2 ms rise
% f_amp = threshold;
% f_tau = round(Fs * .008); %
% f_baselen = round(.005 * Fs);
% % decay
% y = zeros(f_len, 1);
% xxx = [1: f_len];
% f_decaywave = f_amp * exp(xxx./-f_tau);
%
% % rise
% y = zeros(f_rise,1);
% xxx = [1:length(y)];
% f_risewave = (f_amp/length(y))*xxx;
% % plot(f_risewave)
%
% % let's do a 5ms baseline at 0
% f_basewave = zeros(1, f_baselen);
% % f_wave = [f_basewave, f_risewave, f_decaywave];
% f_wave = [f_risewave, f_decaywave];
%
% f_x = [1:length(f_wave)]/Fs;
% plot(f_x,f_wave)
% % fl_fwave = fliplr(f_wave);
% fl_fwave = f_wave;
%
% plot(f_x,fl_fwave)
%
% keyboard
%
% % convolution
% % normI2=normI;
% %maxI2=max(avgwave)*2;
% %normI2(normI2>maxI2)=maxI2;
% conv_PSC = conv(longI_inv, fl_fwave, 'same');
% % conv_sel=conv_PSC(PSC_event);
% % outpoints3=find(conv_sel< -std(conv_sel) );
% % convsel = find(conv_PSC > (2*std(conv_PSC) + mean(conv_PSC) ));
% [pks, PSCidx] = findpeaks(conv_PSC, 'MinPeakHeight', (0.5*std(conv_PSC) + mean(conv_PSC) ));
%
% PSC_event = longT(PSCidx);
% figure
% plot(longT, longI)
% hold on
% scatter(PSC_event, longI(PSCidx), 'o', 'r')

% keyboard
%%
%Now we want to get the windows surrounding the threshold crossings and put
%them in a matrix.
%We will use this matrix to compute template candidates for a filtering
%stage of processing aimed at event detection
%The candidate filters will be generated via either PCA or by taking a mean
%of the detected waveforms


window1=[0.01 0.010]; %First value is for the part before, and the second value is for the values after the threshold crossing
numevents=length(PSC_event_tmp);
window2=round(window1*Fs);
windowlen1=sum(window2) + 1;
eventTvect_tmp = [-window1(1): 1/Fs : window1(2)];

%% THIS IS IMPORTANT - DETERMINING WINDOWS
if EPSCs == 1
    window3=[0.01 0.015]; % THIS IS TO USE AFTER ADJUSTING PEAK - uses a longer window so we can get decay
elseif EPSCs == 0
    window3 = [0.02 0.4];
elseif EPSCs == 2
    window3=[0.01 0.05];
end
window4=round(window3*Fs);
windowlen2= sum(window4) + 1;
eventTvect = [-window3(1): 1/Fs : window3(2)];

% pre-allocate variables
wavemat=nan(numevents,windowlen2);

amplitude(1:numevents) = NaN(1, numevents);
riseT(1:numevents) = NaN(1, numevents);
decayT(1:numevents) = NaN(1, numevents);
tau(1:numevents) = NaN(1, numevents);
taucoeff(1:numevents) = NaN(1, numevents);


win4_1 = window4(1);
win4_2 = window4(2);
win3_1 = window3(1);
% figure; hold on
% keyboard
PSCidxnewnew = PSCidx;
PSC_event_new = PSC_event_tmp;
if slptreadjust == 1
    for ii=1:numevents
        ig=PSCidx(ii);
        d1= ig-window2(1);
        d2= ig+window2(2);
        if d1>0 && d2 < length(normI)
            tmpwav=normI(d1:(d2));
            % now for to grab the amplitude and then its location
            
            tmpamp = max( tmpwav(Fs * .008:Fs*.012) ); % THIS ONLY WORKS IF PRE-WIN IS 10 MS
            
            idx = find(tmpwav == tmpamp);
            
            %% READJUST TIMING AND OFFSET ACTUAL EVENT TIME
            
            PSC_event_new(ii) = PSC_event_tmp(ii) + eventTvect_tmp(idx(1)) ;
            %                     PSCidx_new(ik) = round(eventTvect_tmp(idx(1)) *Fs) + PSCidx(ik);
            PSCidxnewnew(ii) = PSCidx(ii) + round(eventTvect_tmp(idx(1)) *Fs);
            
        end
    end
    
    
    PSC_event = PSC_event_new;
    PSCidx = PSCidxnewnew;
else
    PSC_event = PSC_event_tmp;
end


%%
% grab a templwave wave and scale it real quick from the averages
tmpwavemat=zeros(length(PSC_event),length(eventTvect));
% figure; hold on
for ik=1:length(PSC_event)
    i1=PSC_event(ik)*Fs;
    b1=i1-window4(1);
    b2=i1+window4(2);
    if b1>0 && b2<length(longT)
        tempwav=normI(b1:(b2));
        tmpwavemat(ik,:) = tempwav;
        %    plot(wavemat(ik,:))
    end
    
end
%average and normalize the templatewave
templatewave = mean(tmpwavemat,1);
templatewave = -templatewave/max(templatewave);
templatewave = templatewave - max(templatewave);
templatewave = templatewave/min(templatewave);




%%

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

if poolsize == 0
    poolobj = parpool;
end
%         delete(myCluster.Jobs)

pctRunOnAll warning('off')


% for ik= 1:numevents
%     keyboard
parfor ik=1:numevents
    
    %% NOW REDO ALL THAT SHIT FROM EARLIER IN GRABBING THE WAVE
    i1=PSCidx(ik);
    b1= i1-win4_1;
    b2= i1+win4_2;
    if b1>0 && b2 < length(normI)
        tempwav=normI(b1:b2);
        
        % re-baseline to 0. again. it's not that bad
        tmpI = tempwav(1 : Fs* (win3_1 - .003));
        maxI = max(tmpI);
        minI = min(tmpI);
        histrange = [minI:(maxI - minI)/1000:maxI];
        Ihist = histc(tmpI, histrange);
        [~, idx] = max(Ihist);
        base = histrange(idx);
        %             base = mean(tmpI);
        if ~isempty(base)
            
            tempwav = (tempwav - base);
            
            % recalculate amplitude
            amplitude(ik) = max( tempwav(Fs * (win3_1-.002) :Fs*(win3_1+.002)) ); 
            
            % iN CASE YOU WANT TO PLOT IT
            %            plot(wavemat(ik,:))
            %             plot(eventTvect, tempwav);
            
            
            %%  Now grab rise time, decay, and tau
            idx = round(Fs*win3_1);
            
            % calculate rise time as time for .1 to .9 of peak amplitude
            perc10 = find(tempwav(1:idx(1)) <= .1 * amplitude(ik));
            perc90 = find(tempwav(1:idx(1)) >= .9 * amplitude(ik));
            
            % if we can't get good rise time, don't do it.
            if isempty(perc10) == 1 || isempty(perc90) == 1
                riseT(ik) = NaN;
                decayT(ik) = NaN;
                tau(ik) = NaN;
                taucoeff(ik) = NaN;
            else
                
                riseT(ik) = (perc90(1) - perc10(end))/Fs;
                
                if riseT(ik) < 0
                    riseT(ik) = NaN;
                end
                %                 hold on
                %                 plot(eventTvect(perc90(1)), tempwav(perc90(1)), 'r*')
                %                 plot(eventTvect(perc10(end)), tempwav(perc10(end)), 'r*')
                %                 hold off
                %                 pause
                
                % if we don't get good lambda (decay), then don't do the
                % exponential decay fitting
                Lamp = .37 * amplitude(ik);
                Lampidx = find(tempwav(idx(1):end) <= Lamp);
                
                % get rid of events that are too far off from norm
                % normalize and scale the template wave and the tempwav
                
                normtempwav = tempwav/max( tempwav(Fs * .008:Fs*.012) );
                normtempwav = normtempwav - min(normtempwav);
                normtempwav = normtempwav/max(normtempwav);
                tmpwavesubtracted = normtempwav - templatewave;
                tmpwavesubstd = nanstd(tmpwavesubtracted( 100 :end));
                
                %                         plot(templatewave); hold on; plot(normtempwav); hold off
                %                         pause

                if isempty(Lampidx) == 1 || tmpwavesubstd > 0.15
                    decayT(ik) = NaN;
                    tau(ik) = NaN;
                    taucoeff(ik) = NaN;
                else
                    decayT(ik) = Lampidx(1)/Fs;
                    
                    % then calculate tau!!!!! like real decay bro
                    
                    tofity = tempwav(idx(1):end)';
                    tofity = tofity/tofity(1);
                    tofitx = eventTvect(idx(1):end)';
                    tofitx = tofitx - tofitx(1);
                    
                    % fitting double exponential f2(x) = a*exp(b*x) + c*exp(d*x)
                    % but here let's just do a single exponential
                    %                             ft = fittype( 'a*exp(-x*b)', 'independent', 'x', 'dependent', 'y' );
                    %                             taufit = fit(tofitx , tofity, ft, 'StartPoint', [tofity(1), decayT(ik)],'MaxIter',1000);
                    try
                        taufit = fit(tofitx , tofity, 'exp2');
                        % for weighted tau: (A1*t1 + A2*t2)/ (A1 +A2)
                        % from  f2(x) = a*exp(b*x) + c*exp(d*x)
                        tau(ik) = ( (taufit.a*(1/taufit.b)) + (taufit.c*(1/taufit.d)) ) / (taufit.a + taufit.c)  ;
                        % t1*a1+ t2 *(1-a1)
                        %                             tau(ik) = (-taufit.a*taufit.b) + (taufit.d*(1-taufit.a))  ;
                        tau(ik) = abs( tau(ik) ) ;
                        taucoeff(ik) = (taufit.a + taufit.c);
                        
                        if tau(ik) > 0.3
                            tau(ik) = NaN;
                            taucoeff(ik) = NaN;
                        end
                        
                        % store the wave into 'wave'
                        wavemat(ik,:) = tempwav;
                        
                    catch
                        tau(ik) = NaN;
                        taucoeff(ik) = NaN;
                    end
                    %
                    %                                                         plot(taufit,tofitx,tofity)
                    %                                                         title(['tau = ' num2str(tmptau)])
                    %                                                         pause
                    %                                                         hold on
                    %                             %                             scatter(decayT(ik), Lamp, 'o', 'g')
                    %                                                         hold off
                    %                             keyboard
                end
            end
        end
    end
end



dataPSC.amplitude = amplitude;
dataPSC.riseT = riseT;
dataPSC.decayT = decayT;
dataPSC.tau = tau;
dataPSC.taucoeff = taucoeff;

end
