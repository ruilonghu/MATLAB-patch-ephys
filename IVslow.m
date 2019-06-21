% VC IV
warning('off')
close all
clearvars -except Trace_*

savename =  '190619 slice1 (1) AOB GC' ;
yes_save =   1  ;

% which experiment
Experimentnum =            1            ;
% load which trace?
Trace =                    3            ;
% channel traces for current and voltage
curr_channel =  1   ;
volt_channel =  2   ;

%% params
baseline_duration = .5; % seconds
stimulus_duration = 1; % seconds
stimulus_end = baseline_duration + stimulus_duration;

Experiment = ['Trace_' num2str(Experimentnum) '_'];
tracelist_curr = who([Experiment, num2str(Trace),'_*_', num2str(curr_channel)]);
tracelist_volt = who([Experiment, num2str(Trace),'_*_', num2str(volt_channel)]);

% reorders tracelist if it goes above 10
for i = 1:length(tracelist_curr)
    reord_tracelist_current{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(curr_channel)];
end

for i = 1:length(tracelist_volt)
    reord_tracelist_voltage{i} = [Experiment, num2str(Trace), '_', num2str(i), '_', num2str(volt_channel)];
end

% gets sampling rate
tmp = eval(reord_tracelist_voltage{1});
Fs = (length(tmp(:,1))/max(tmp(:,1)));

for i = 1:length(tracelist_curr)
    
    [volt_trace] = eval(reord_tracelist_voltage{i});
    [curr_trace] = eval(reord_tracelist_current{i});
    
    data_T(:,i) = curr_trace(:,1);
    data_I(:,i) = curr_trace(:,2);
    data_V(:,i) = volt_trace(:,2);
    
    % uses mode instead of mean
    base_curr_mean = mode(curr_trace((1:((baseline_duration - .2)*Fs)), 2));
    base_volt_mean = mode(volt_trace((1:((baseline_duration - .2)*Fs)), 2));
    % voltage step
    idxstart = ((baseline_duration * Fs) + 1);
    idxend = ((stimulus_end * Fs) - 1);
    voltages(i) = mean(volt_trace(idxstart:idxend, 2));
    voltage_steps(i) = voltages(i)- base_volt_mean;
    
    % current step
    idxstart = ((stimulus_end - .02)* Fs);
    idxend = (stimulus_end * Fs);
    current_steps(i) = mean(curr_trace(idxstart:idxend, 2));
    
    %     subtracted_curr_steps(i) = current_steps(i) - (base_currstep* i* (10^12));
    subtracted_curr_steps(i) = current_steps(i)- (base_curr_mean);
    subtracted_curr_pA(i) = subtracted_curr_steps(i)  * 10^12;
    
    voltage_steps_mV(i) = (voltage_steps(i) * 10^3);
    
    R_MOhm(i) = abs(voltage_steps(i)/subtracted_curr_steps(i) * 10^-6); % converts to mega ohms
    voltages_mV(i) =  round((voltages(i) * 10^3));
end

%% plot IV

IVfig=figure;

subplot(2,1,1)

plot(data_T, data_I)

subplot(2,1,2)
plot(voltages_mV, subtracted_curr_pA, 'o-')
xlabel('voltage (mV)', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold')
ylabel('current (pA)','FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

box off;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 10 7.5]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 7.5]);
set(gcf, 'PaperPosition', [0 0 10 7.5]);




%% save data

if yes_save == 1
    saveas(IVfig, [savename ' IVslow.jpg'], 'jpg')
    save([savename, ' IVslow.mat'], '-regexp',  '^(?!Trace_.*$).')
end