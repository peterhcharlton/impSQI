function res = imp_rr_sqi(sig, fs)
% imp_rr_sqi estimates respiratory rates (RRs) from impedance pneumography
% signals and calculates a corresponding signal quality index (SQI).
%
%               imp_rr_sqi(sig, fs)
%
%	Inputs:
%       sig - an impedance pneumography signal in a vector of sample values.
%       fs - the sampling frequency of the signal (in Hz)
%
%	Outputs:
%       res - a structure containing the following fields (with one value per window):
%               - flat_line: whether or not this window was a flat line
%               - sqi_novel: the SQI value (1 = high quality, 0 = low quality)
%               - rr_novel: the estimated respiratory rate (RR) using the technique of the novel SQI
%               - prop_norm_dur: proportion of window taken up by normal breath durations
%               - prop_bad_breaths: proportion of breaths which have outlying breath durations
%               - R2: mean correlation coefficient
%               - R2min: normalised standard deviation of breath durations
%               - win_no: window number
%               - win_start_sample: first sample of this window
%               - win_end_sample: last sample of this window
%           
%   Further Information:
%       This version of imp_rr_sqi is intended to be a reproduction of the algorithm
%       reported in:
%           Charlton P. H. et al., "An impedance pneumography signal
%           quality index for respiratory rate monitoring: design,
%           assessment and application", https://doi.org/10.1016/j.bspc.2020.102339
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/imp_sqi.html
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.0.1 - 24th Oct 2024 by Peter H Charlton
%
%   Licence:
%       This program is available under the GNU public license. This
%       program is free software: you can redistribute it and/or modify 
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version. This program is distributed in
%       the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%       even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%       PARTICULAR PURPOSE.  See the GNU General Public License for more
%       details: <http://www.gnu.org/licenses/>.
%

%% Setup
% Declare Universal Parameters for use in the algorithm
up = setup_universal_params;
% Re-format variables
imp = reformat_vars(sig,fs);

%% Detecting breaths
res = calc_rr_and_sqi(imp,up);

end

function up = setup_universal_params

% whether to use additional methods for RR estimation and SQI calculation
up.paramSet.do_additional_methods = false;

% window parameters
up.paramSet.winLeng = 32;     % the duration of each window in secs
up.paramSet.winStep = 0; %up.paramSet.winLeng/2; %3;      % the number of secs between each consecutive window

% Eliminate HFs (above resp freqs)
up.paramSet.elim_hf.Fpass = 1.2;  % in Hz
up.paramSet.elim_hf.Fstop = 0.895;  % in Hz
up.paramSet.elim_hf.Dpass = 0.057501127785;
up.paramSet.elim_hf.Dstop = 0.01;

% duration of Tukey window taper in secs
up.paramSet.tukey_win_duration_taper = 2;

% Range of plausible RR freqs
up.paramSet.rr_range = [4 60];

% plotting
up.settings.do_plot = false;
% save folder (same as directory of this script)
mfilePath = mfilename('fullpath');
[up.paths.plots_save_folder,name,ext] = fileparts(mfilePath);

end

function imp = reformat_vars(sig,fs)

% Re-format the variables in the same way as in the original code

imp.v = sig(:);
imp.fs = fs;

end

function res = calc_rr_and_sqi(imp,up)

%% Extract time vector
imp.t = [0:(length(imp.v)-1)]/imp.fs;

%% Filter impedance signal
% Filtering impedance signal to eliminate high frequencies above respiration
imp_filt = lpf_to_exclude_resp(imp, up);

%% Make plot of the impedance signal
if up.settings.do_plot

    figure('Position', [20,20,800,300])
    ftsize = 20; lwidth = 2;
    plot(imp.t, imp.v),
    ylab = ylabel('ImP', 'FontSize', ftsize, 'Rotation', 0);
    xlabel('Time [s]', 'FontSize', ftsize)
    set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    set(gca, 'FontSize', ftsize)
    xlim([0 10])
    
    savepath = [up.paths.plots_save_folder, filesep, 'sample_recording'];
    PrintFigs(gcf, savepath)
    close all
    
end

%% Normalise signals in each window
% The impedance signal corresponding to each window is extracted, then
% it is normalised to have a mean of 0 and standard deviation of 1.
win_starts = imp.t(1):up.paramSet.winLeng:imp.t(end);
[res.flat_line, res.sqi_novel, res.rr_novel, res.prop_norm_dur, res.prop_bad_breaths, res.R2, res.R2min, res.win_no, res.win_start_sample, res.win_end_sample] = deal(nan(length(win_starts),1));
if up.paramSet.do_additional_methods
    [res.rr_agree, res.rr_wch, res.rr_cto] = deal(nan(length(win_starts),1));
end
for win_no = 1 : length(win_starts)

    % input subject
    res.win_no(win_no) = win_no;

    % select data for this window
    win_start = win_starts(win_no);
    win_end = win_start + up.paramSet.winLeng;
    rel_els = find(imp_filt.t >= win_start & imp_filt.t <= win_end);
    rel_data.t = imp_filt.t(rel_els);
    rel_data.v = imp_filt.v(rel_els);

    % store start and end els of this window
    res.win_start_sample(win_no) = rel_els(1);
    res.win_end_sample(win_no) = rel_els(end);

    % determine whether it was a flat line - if so then skip this window
    if length(unique(imp.v(rel_els))) == 1 || sum(isnan(imp.v(rel_els)))>=1
        res.flat_line(win_no) = 1;
        continue
    end
    res.flat_line(win_no) = 0;
    % interpolate data to be at fixed time values
    downsample_freq = 5;
    interp_data.t = win_start: (1/downsample_freq) : win_end;
    interp_data.v = interp1(rel_data.t, rel_data.v, interp_data.t, 'linear');
    % get rid of nans
    interp_data.v(isnan(interp_data.v)) = median(interp_data.v(~isnan(interp_data.v)));
    % normalise data
    rel_imp.v_n = (interp_data.v - mean(interp_data.v))/std(interp_data.v);
    rel_imp.t_n = interp_data.t;
    rel_imp.fs = downsample_freq;
    clear interp_data rel_data downsample_freq rel_els

    %% Estimate RR from the impedance signal

    % - using the novel modified count-orig method (which also provides signal quality metrics)
    [res.sqi_novel(win_no), res.rr_novel(win_no), res.prop_norm_dur(win_no), res.prop_bad_breaths(win_no), res.R2(win_no), res.R2min(win_no)] ...
        = ref_cto_mod(rel_imp, up, 'no');
    
    if up.paramSet.do_additional_methods

        % - using the original count-orig method (returns a nan if poor-quality)
        res.rr_cto(win_no,1) = ref_cto(rel_imp, up);

        % - using an FFT
        res.rr_wch(win_no,1) = ref_wch(rel_imp, up);

        % using the agreement SQI
        res.rr_agree(win_no,1) = mean([res.rr_cto(win_no), res.rr_wch(win_no)]);

    end

    clear rel_els rel_data win_start win_end rel_imp

end
clear imp imp_filt no_wins qual win_starts win_no


%% Calculate agreement SQI
if up.paramSet.do_additional_methods

    good_els = abs(res.rr_cto-res.rr_wch) < 2;
    res.agree_sqi = false(length(res.rr_cto),1);
    res.agree_sqi(good_els) = true;
    clear good_els

end

end

function imp = lpf_to_exclude_resp(imp, up)

imp.v(isnan(imp.v)) = median(imp.v(~isnan(imp.v)));

%% Window signal to reduce edge effects
duration_of_signal = imp.t(end) - imp.t(1);
prop_of_win_in_outer_regions = 2*up.paramSet.tukey_win_duration_taper/duration_of_signal;
tukey_win = tukeywin(length(imp.v), prop_of_win_in_outer_regions);
d_s_win = imp;    % copy time and fs
d_s_win.v = detrend(imp.v(:)).*tukey_win(:);

%% LPF to remove freqs above resp
respWave.t = d_s_win.t;
respWave.v = lp_filter_signal_to_remove_freqs_above_resp(d_s_win.v, d_s_win.fs, up);
respWave.fs = d_s_win.fs;
imp = respWave;
end

function s_filt = lp_filter_signal_to_remove_freqs_above_resp(s, Fs, up)
%% Filter pre-processed signal to remove freqs above resp

% parameters for the low-pass filter to be used
flag  = 'scale';
Dpass = up.paramSet.elim_hf.Dpass;
Dstop = up.paramSet.elim_hf.Dstop;
Fstop = up.paramSet.elim_hf.Fstop;
Fpass = up.paramSet.elim_hf.Fpass;

% create filter
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [1 0], [Dstop Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0159;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(Fs/2);

% Prepare signal
s_dt=detrend(s);
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
end

function [qual, rr_cto, prop_norm_dur, prop_bad_breaths, R2, R2min] = ref_cto_mod(sum_both, up, save_name)

sum_both.t_n = sum_both.t_n-sum_both.t_n(1);
sum_both.v_n = -1*detrend(sum_both.v_n);

%% Identify relevant peaks and troughs

% identify peaks
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% identify troughs
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% define peaks threshold
q3 = quantile(sum_both.v_n(peaks), 0.75);
thresh = 0.2*q3;
% find relevant peaks
rel_peaks = peaks(sum_both.v_n(peaks) > thresh);
% define troughs threshold
q3t = quantile(sum_both.v_n(troughs), 0.25);
thresh = 0.2*q3t;
% find relevant troughs
rel_troughs = troughs(sum_both.v_n(troughs) < thresh);

%% find valid breathing cycles
% exclude peaks which aren't the highest between a pair of consecutive
% troughs:
invalid_peaks = zeros(length(rel_peaks),1);
for trough_pair_no = 1 : (length(rel_troughs)-1)
    
    % identify peaks between this pair of troughs
    cycle_rel_peak_els = find(rel_peaks > rel_troughs(trough_pair_no) & rel_peaks < rel_troughs(trough_pair_no+1));
    cycle_rel_peaks = rel_peaks(cycle_rel_peak_els);
    if length(cycle_rel_peaks) > 1
        [~, rel_el] = max(sum_both.v_n(cycle_rel_peaks));
        bad_rel_peaks_els = setxor(1:length(cycle_rel_peak_els), rel_el);
        invalid_peaks(cycle_rel_peak_els(bad_rel_peaks_els)) = 1;
    end
end
rel_peaks = rel_peaks(~invalid_peaks);

% if there is more than one initial peak (i.e. before the first trough) then take the highest:
initial_peaks = find(rel_peaks < rel_troughs(1));
other_peaks = find(rel_peaks >= rel_troughs(1));
if length(initial_peaks)>1
    [~, rel_initial_peak] = max(sum_both.v_n(rel_peaks(initial_peaks)));
    rel_peaks = rel_peaks([rel_initial_peak, other_peaks]);
end

% valid cycles start with a peak:
valid_cycles = false(length(rel_peaks)-1,1);
cycle_durations = nan(length(rel_peaks)-1,1);

for peak_no = 2 : length(rel_peaks)
    
    % exclude if there isn't a rel trough between this peak and the
    % previous one
    cycle_rel_troughs = rel_troughs(rel_troughs > rel_peaks(peak_no-1) & rel_troughs < rel_peaks(peak_no));
    if length(cycle_rel_troughs) ~= 0
        valid_cycles(peak_no-1) = true;
        cycle_durations(peak_no-1) = sum_both.t_n(rel_peaks(peak_no)) - sum_both.t_n(rel_peaks(peak_no-1));
    end
end
valid_cycle_durations = cycle_durations(valid_cycles);

% Calc RR
if isempty(valid_cycle_durations)
    rr_cto = nan;
else
    % Using average breath length
    ave_breath_duration = mean(valid_cycle_durations);
    rr_cto = 60/ave_breath_duration;
end

%% Resiratory SQI

if isnan(rr_cto)
    qual = false;
    prop_norm_dur = 0;
    prop_bad_breaths = 100;
    R2 = 0;
    R2min = 0;
else
    
    %find mean breath-to-breath interval to define size of template
    rr=floor(mean(diff(rel_peaks)));
    ts=[];
    j=find(rel_peaks>rr/2);
    l=find(rel_peaks+floor(rr/2)<length(sum_both.v_n));
    new_rel_peaks = rel_peaks(j(1):l(end));
    if isempty(new_rel_peaks)
        qual = false;
        prop_norm_dur = 0;
        prop_bad_breaths = 100;
        R2 = 0;
        R2min = 0;
        return
    else
        %find breaths
        for k=1:length(new_rel_peaks)
            t=sum_both.v_n(new_rel_peaks(k)-floor(rr/2):new_rel_peaks(k)+floor(rr/2));
            tt=t/norm(t); tt = tt(:)';
            ts=[ts;tt];
        end
    end
    
    %find ave template
    if size(ts,1) > 1
        avtempl=mean(ts,1);
    else
        avtempl=nan(size(ts));
    end
    
    %now calculate correlation for every beat in this window
    r2 = nan(size(ts,1),1);
    for k=1:size(ts,1)
        r2(k)=corr2(avtempl,ts(k,:));
    end
    %calculate mean correlation coefficient
    R2=mean(r2);
    R2min = std(valid_cycle_durations)/mean(valid_cycle_durations);
    %peak_heights = sum_both.v_n(rel_troughs);
    %R2min = std(peak_heights)/mean(peak_heights);
    
    % calculate number of abnormal breath durations
    median_dur = median(valid_cycle_durations);
    temp = valid_cycle_durations > (1.5*median_dur) | valid_cycle_durations < (0.5*median_dur);
    prop_bad_breaths = 100*sum(temp)/length(temp);
    
    % find prop of window taken up by normal breath durations
    norm_dur = sum(valid_cycle_durations(~temp));
    win_length = sum_both.t_n(end) - sum_both.t_n(1);
    prop_norm_dur = 100*norm_dur/win_length;
    
    % determine whether this window is high or low quality
    if prop_norm_dur > 60 && prop_bad_breaths < 15 && R2 >= 0.75 && R2min < 0.25
        qual = true;
    else
        qual = false;
    end
    
end

%% Plot template and inidividual beats
save_name = 'sqi_process';
save_name = 'no';
if ~strcmp(save_name, 'no') && ~isnan(rr_cto)
    
    paper_size = [12, 8];
    figure('Position', [50, 50, 100*paper_size(1), 100*paper_size(2)], 'Color',[1 1 1])
    lwidth1 = 3; lwidth2 = 2; ftsize = 18;
    % plot signal
    subplot(2,2,[1,2]), plot(sum_both.t_n-sum_both.t_n(1), sum_both.v_n, 'LineWidth', lwidth2), hold on
    plot(sum_both.t_n(rel_peaks(logical([valid_cycles; 1])))-sum_both.t_n(1), sum_both.v_n(rel_peaks(logical([valid_cycles; 1]))), '.r', 'MarkerSize', 20)
    %plot(sum_both.t_n(new_rel_peaks)-sum_both.t_n(1), sum_both.v_n(new_rel_peaks), '.r', 'MarkerSize', 20)
    %plot(sum_both.t_n(rel_troughs)-sum_both.t_n(1), sum_both.v_n(rel_troughs), '.k', 'MarkerSize', 20)
    xlim([0, sum_both.t_n(end)-sum_both.t_n(1)])
    xlabel('Time [s]', 'FontSize', ftsize)
    ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    set(gca, 'FontSize', ftsize, 'YTick', [])
    % plot template
    time = 0:(length(avtempl)-1); time = time./sum_both.fs;
    subplot(2,2,3), hold on,
    for beat_no = 1 : size(ts,1)
        plot(time, ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
    end
    plot(time, avtempl, 'r', 'LineWidth', lwidth1)
    set(gca, 'YTick', [])
    xlabel('Time [s]', 'FontSize', ftsize)
    xlim([0, time(end)])
    ylab=ylabel('ImP', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    set(gca, 'FontSize', ftsize)
    %set ylim
    rang = range(ts(:));
    ylim([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]);
    % prop_norm_dur > 60 && prop_bad_breaths < 15 && R2 >= 0.75 && R2min < 0.25
    if R2 >= 0.75, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.25, 0.1,0.1],'String',['R^2 = ' num2str(R2, 2)], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    if prop_bad_breaths < 15, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.18, 0.1,0.1],'String',['prop invalid breaths = ' num2str(prop_bad_breaths, 2), ' %'], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    if prop_norm_dur > 60, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.11, 0.1,0.1],'String',['prop valid duration = ' num2str(prop_norm_dur, 2), ' %'], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    if R2min < 0.25, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.04, 0.1,0.1],'String',['Breath interval variability = ' num2str(100*R2min, 2), ' %'], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    
    %annotation('textbox',[0.5, 0.1, 0.1,0.1],'String',{['R^2 = ' num2str(R2, 2)] , ['prop breaths bad = ' num2str(prop_bad_breaths, 2) '%'], ['prop dur good = ' num2str(prop_norm_dur,2) '%'], ['norm SD durations = ' num2str(R2min,2) '%']}, 'FontSize', ftsize, 'LineStyle', 'None')
    if qual
        annotation('textbox',[0.8, 0.15, 0.1,0.1],'String','High Quality', 'Color', 'b', 'FontSize', ftsize+4, 'LineStyle', 'None')
    else
        annotation('textbox',[0.8, 0.15, 0.1,0.1],'String','Low Quality', 'Color', 'r', 'FontSize', ftsize+4, 'LineStyle', 'None')
    end
    savepath = [up.paths.plots_save_folder, save_name];
    PrintFigs(gcf, paper_size, savepath, up)
    close all
end

%% Plot summary figure
save_name = 'resp_sqi_summary_viva';
save_name = 'no';
if ~strcmp(save_name, 'no') && ~isnan(rr_cto)
    
    % check to see if a figure has been opened
    paper_size = [16, 7];
    if isempty(findall(0,'Type','Figure'))
        figure('Position', [50, 50, 100*paper_size(1), 100*paper_size(2)], 'Color',[1 1 1])
    end
    lwidth1 = 3; lwidth2 = 2; ftsize = 18; plotted = 0;
    
    % see if this signal is appropriate
    if R2 < 0.55
        % Plot signal
        subplot(2,3,1:2),plot(sum_both.t_n-sum_both.t_n(1), -1*sum_both.v_n, 'LineWidth', lwidth2), hold on
        plot(sum_both.t_n(rel_peaks(logical([valid_cycles; 1])))-sum_both.t_n(1), -1*sum_both.v_n(rel_peaks(logical([valid_cycles; 1]))), '.k', 'MarkerSize', 20), hold off
        plotted = 1;
        xlim([0, sum_both.t_n(end)-sum_both.t_n(1)]), ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); set(gca, 'FontSize', ftsize, 'YTick', [])
        title('Low quality signal segment', 'FontSize', ftsize)
        %set ylim
        rang = range(sum_both.v_n);
        ylim(-1*fliplr([min(sum_both.v_n(:))-0.1*rang, max(sum_both.v_n(:))+0.1*rang]));
        % plot template
        subplot(2,3,3), cla, hold on
        time = 0:(length(avtempl)-1); time = time./sum_both.fs;
        for beat_no = 1 : size(ts,1)
            plot(time, -1*ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
        end
        plot(time, -1*avtempl, 'r', 'LineWidth', lwidth1), hold off
        set(gca, 'YTick', [])
        xlim([0, 6]) % Jan Submission
        xlim([0, time(end)])
        title(['Template (R^2 = ' num2str(R2,2) ')'], 'FontSize', ftsize)
        set(gca, 'FontSize', ftsize)
        %set ylim
        rang = range(ts(:));
        ylim(-1*fliplr([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]));
    elseif R2 >= 0.97
        % Plot signal
        subplot(2,3,4:5),plot(sum_both.t_n-sum_both.t_n(1), -1*sum_both.v_n, 'LineWidth', lwidth2), hold on
        plot(sum_both.t_n(rel_peaks(logical([valid_cycles; 1])))-sum_both.t_n(1), -1*sum_both.v_n(rel_peaks(logical([valid_cycles; 1]))), '.k', 'MarkerSize', 20), hold off
        plotted = 1;
        xlim([0, sum_both.t_n(end)-sum_both.t_n(1)]), xlabel('Time [s]', 'FontSize', ftsize), ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); set(gca, 'FontSize', ftsize, 'YTick', [])
        title('High quality signal segment', 'FontSize', ftsize)
        %set ylim
        rang = range(-1*sum_both.v_n);
        ylim(-1*fliplr([min(-1*sum_both.v_n(:))-0.1*rang, max(-1*sum_both.v_n(:))+0.1*rang]));
        
        % plot template
        subplot(2,3,6), cla, hold on
        time = 0:(length(avtempl)-1); time = time./sum_both.fs;
        for beat_no = 1 : size(ts,1)
            plot(time, -1*ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
        end
        plot(time, -1*avtempl, 'r', 'LineWidth', lwidth1), hold off
        set(gca, 'YTick', [])
        xlabel('Time [s]', 'FontSize', ftsize)
        xlim([0, 3]) % Jan Submission
        xlim([0, time(end)])
        %ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
        %set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        title(['Template (R^2 = ' num2str(R2,2) ')'], 'FontSize', ftsize)
        set(gca, 'FontSize', ftsize)
        %set ylim
        rang = range(ts(:));
        ylim(-1*fliplr([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]));
    end
        
        
%     annotation('textbox',[0.5, 0.1, 0.1,0.1],'String',{['R2 = ' num2str(R2, 2)] , ['prop breaths bad = ' num2str(prop_bad_breaths, 2) '%'], ['prop dur good = ' num2str(prop_norm_dur,2) '%'], ['norm SD durations = ' num2str(R2min,2) '%']}, 'FontSize', ftsize, 'LineStyle', 'None')
%     if qual
%         annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','High Quality', 'Color', 'b', 'FontSize', ftsize, 'LineStyle', 'None')
%     else
%         annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','Low Quality', 'Color', 'r', 'FontSize', ftsize, 'LineStyle', 'None')
%     end
    savepath = [up.paths.plots_save_folder, save_name];
    %PrintFigs(gcf, paper_size, savepath, up)
    %close all
end


end

function rr_cto = ref_cto(sum_both, up)

% identify peaks
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% identify troughs
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% define threshold
q3 = quantile(sum_both.v_n(peaks), 0.75);
thresh = 0.2*q3;
% find relevant peaks and troughs
extrema = sort([peaks(:); troughs(:)]);
rel_peaks = peaks(sum_both.v_n(peaks) > thresh);
rel_troughs = troughs(sum_both.v_n(troughs) < 0);

% find valid breathing cycles
% valid cycles start with a peak:
valid_cycles = zeros(length(rel_peaks)-1,1);
cycle_durations = nan(length(rel_peaks)-1,1);
for peak_no = 1 : (length(rel_peaks)-1)
    
    % valid if there is only one rel_trough between this peak and the
    % next
    cycle_rel_troughs = rel_troughs(rel_troughs > rel_peaks(peak_no) & rel_troughs < rel_peaks(peak_no+1));
    if length(cycle_rel_troughs) == 1
        valid_cycles(peak_no) = 1;
        cycle_durations(peak_no) = sum_both.t_n(rel_peaks(peak_no+1)) - sum_both.t_n(rel_peaks(peak_no));
    end
    
end

% Calc RR
if sum(valid_cycles) == 0
    rr_cto = nan;
else
    % Using average breath length
    ave_breath_duration = nanmean(cycle_durations);
    rr_cto = 60/ave_breath_duration;
end

end

function rr_wch = ref_wch(sum_both, up)

segLen = 2^nextpow2(12*sum_both.fs);
noverlap = segLen/2;
[data.power, data.freqs] = pwelch(sum_both.v_n,segLen,noverlap, [], sum_both.fs);

% Find spectral peak
[rr_wch, ~, ~] = find_spectral_peak(data, up);

end

function [v, f, p] = find_spectral_peak(data, up)
freq_range = up.paramSet.rr_range/60;

cand_els = zeros(length(data.power),1);
for s = 2 : (length(data.power)-1)
    if data.power(s) > data.power(s-1) & data.power(s) > data.power(s+1) & data.freqs(s) > freq_range(1) & data.freqs(s) < freq_range(2)
        cand_els(s) = 1;
    end
end
clear s freq_range

temp = data.power - min(data.power);  % because the FFT and ACF have -ve values for the power spectra.
[~, r_el] = max(temp.*cand_els); clear cand_els
r_freq = data.freqs(r_el); clear r_el
v = 60*r_freq;

%% Store spectrum
f = data.freqs;
p = data.power;
end

function PrintFigs(h, savepath, do_close)

fprintf('\n - Saving: %s', savepath)

print(h,'-dpng',savepath)

if ~(nargin == 3 && do_close == 0)
    close all
end

end
