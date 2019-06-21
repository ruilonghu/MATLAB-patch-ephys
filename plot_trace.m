% function plottrace = plot_trace(trace_list)
close all
clearvars -except Trace_*

% plot traces
% which experiment

% which experiment 
Experimentnum =              1              ;
% load which  traces?      
Trace =    [                 8            ]; %15 16 17
   
Experiment = ['Trace_' num2str(Experimentnum) '_'];

% set channel
channel =           1              ; 

% offset? 
offset =            0              ;

do_axis =              0              ;
yes_save =             0              ;
TITLE = '190517 slice8 (17) MOB FS CsClQX - MCPO-ChR2 2ms @1s';
% MCPO-ChR2 2ms @1s
% Thy1ChR2 -70mV 2ms
% IVfast -60mV
saveascii = 0; 


%%

figure; hold on
if offset == 1
    
    tmp = eval([Experiment, num2str(Trace(1)), '_1_', num2str(1)]);
    maxT = round(tmp(end,1));
    Fs = length(tmp(:,1))/maxT;
    
    data_T = [];
    data_V = [];
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
        data_V = vertcat(data_V, longV);
        
        lastTi = data_T(end);
    end
%%

%     [bLP,aLP] = butter(2,[5 200]/(Fs/2),'bandpass'); %300Hz LP filter 2nd order
% 
% %     [bLP,aLP] = butter(2, 200/(Fs/2), 'low'); %300Hz LP filter 2nd order
%     data_V = (filtfilt(bLP,aLP, double(data_V)));
%     
    % figure
    plot(data_T, data_V, 'k')
    
    %         set(gca, 'FontSize', 14, 'FontName', 'Arial','FontWeight', 'bold');
    
    if do_axis == 1
        
        box off;
        hold on;
        set(gcf,'Color',[1 1 1]);
        set(gcf,'Units','inches');
        set(gcf,'Position',[1 1 10 7.5]);
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [10 7.5]);
        set(gcf, 'PaperPosition', [0 0 10 7.5]);
        axis off
    end

    %% plotting without offset
elseif offset == 0
    
    prevlen = 0;
    
    for i = 1:length(Trace)
        tracelist = who([Experiment, num2str(Trace(i)), '_*_', num2str(channel)]);
       for j = 1:length(tracelist)
                reord_tracelist{j} = [Experiment,num2str(Trace(i)), '_', num2str(j), '_', num2str(channel)];
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
            
            data_V(:, ( ((i-1) * prevlen) +j)) = trace(:,2);

            plot(trace(:,1), trace(:,2));
            
            
        end
        prevlen = length(tracelist);
    end
    
end
%         set(gca, 'FontSize', 18, 'FontName', 'Arial','FontWeight', 'bold');


% stdataV = std(data_V, 1);


if do_axis == 1
    
    box off;
    hold on;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','inches');
    set(gcf,'Position',[1 1 10 7.5]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 7.5]);
    set(gcf, 'PaperPosition', [0 0 10 7.5]);
    
    axis off
end


if yes_save == 1
    plotplot = gca;
    saveas(plotplot, [TITLE '.jpg'])
    save([TITLE, ' plottrace.mat'], '-regexp',  '^(?!Trace_.*$).')
    
    
%     print -depsc -tiff -r300 -painters STA2.eps
end



if saveascii == 1
    asciisavedat = horzcat(data_T,data_V);
    csvwrite([TITLE,'_dat.csv'], asciisavedat)
%     csvwrite([TITLE,'_dat.csv'], data_V)
end