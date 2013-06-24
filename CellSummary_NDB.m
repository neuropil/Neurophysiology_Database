%% CellSummary_gt
% Make a raster with solo and nspike data.
%
% USAGE: [fig_handle, raster_info] = CellSummary_gt(behav_file, spk_file, align_ind, window,...
%      by_mixture_ratio_flag, trials_to_include)
% EXAMPLE: h = CellSummary_gt(behav_fname, spk_fname, 3, [-2 1], 0, [1:145]);
%
% INPUTS:  behav_file - path and filename of ratbase file
%          spk_file - path and filename of spike file
%          align_ind - which trial event to align rasters/psth to (indicated in 
%              GetSCGlobals.m)
%               OdorPokeIn = 1
%               OdorValveOn = 2
%               OdorPokeOut = 3
%               WaterPokeIn = 4
%               WaterValveOn = 5
%               WaterPokeOut = 6
%               NextOdorPokeIn = 7
%          window - 2x1 vector of times to rasterize spikes, relative to
%              event_times ([START STOP]) (e.g., [-1 1])
%          by_mixture_ratio_flag - whether to consider separate mixture
%              ratios to be separate stimuli, or not (1 - yes; 0 - no)
%          trials_to_include - which behavioral trials to include in
%              analysis (OPTIONAL, default includes all trials)
%
% OUTPUTS: fig_handle - handle to figure of rasters, psths
%          raster_info - structure containing info capable of
%              reconstructing rasters
%
% Adapted from CellSummary for use w/ go tone data.  Differs in that trials are grouped by
% response (L/R), mixture ratio, and go tone
% delay.
%
% 5/19/08 - GF
% Updates:
% $ 12/24/08 - added by_mixture_ratio flag - can be set to 0 for raster
%      display, 1 for calculation of raster_info

function [raster_info] = CellSummary_NDB(behav_file, spk_file, align_ind, window,...
    by_mixture_ratio_flag, trials_to_include)

% MIN_TRIAL_PLOTS = 6; % for psths
%PSTH_SMOOTH_FACTOR = (window(2) - window(1)) * 10; % this is the "standard value" as of 8/2011
PSTH_SMOOTH_FACTOR = (window(2) - window(1)) * 5;

%% load spike file

spike_times = spk_file/1000000; % converts to time in seconds

odorid = [95 5 60 40 80 20 50];

%% load behavioral file

before_trial_res = length(behav_file.start);

if exist('trials_to_include','var') % don't include all trials
    behav_file = RestrictSession(behav_file, trials_to_include);
end

after_trial_res = length(behav_file.start);

GetNeuroGlobals;
global RESOLUTION; % set to 1000 to plot spike times to 1 msec resolution
global TRIAL_EVENTS;
trial_events = fieldnames(TRIAL_EVENTS);
trial_events{2} = 'DIO'; %DIO is OdorValveOn
trial_events{5} = 'WaterDeliv'; %WaterDeliv is WaterValveOn
trial_events{8} = 'cue'; %WaterDeliv is WaterValveOn


% PresentationFigSetUp; %%% THIS IS WHERE trial_events comes from (weird)
% fig_handle = gcf;




%% get event timings of each trial type

if by_mixture_ratio_flag == 0
       [raster_events, stim_params] = TrialClassification_dispatcher2test(behav_file, trial_events);
%     temp = 0;
elseif by_mixture_ratio_flag == 1
       [raster_events, stim_params] = TrialClassification_dispatcher_previousTrial(behav_file, trial_events);
%     temp = 1;
else
    error('by_mixture_ratio_flag must be 0 or 1');
end

%% determine the total number of (performed) trials
num_stim = length(raster_events.L.trial_start);
total_num_trials = 0; % initialize
for stim_num = 1:num_stim
    total_num_trials = total_num_trials + length(raster_events.L.trial_start{stim_num}) + length(raster_events.R.trial_start{stim_num});
end

%% Raster setup constants
raster_l = 0.08;
raster_w = 0.8;
raster_u_buff = 0.033;
total_rasters_h = 0.72; % fraction of total figure height to devote to rasters
intrastim_space = 0.005;
intrachoice_space = intrastim_space / 3;
available_rasters_h = total_rasters_h - (((num_stim - 1) * intrastim_space) + (num_stim * intrachoice_space)); 
textbuff = 0.04;


%% choice bar setup constants
choicebar_hor_buff = 0.01;
choicebar_w = 0.005;
choicebar_l = raster_l - (choicebar_hor_buff + choicebar_w);

%% stimulus bar setup constants
% stimbar_hor_buff = choicebar_hor_buff;
% stimbar_w = choicebar_w;
% stimbar_l = raster_l + raster_w + stimbar_hor_buff;

raster_u = 1 - raster_u_buff; % initialize

SECONDARY_EVENTS_COLORS = [...
    0 1 0;...         %  Neon Green OdorPokeIn
    1 0.5 1;...       %  Hot Pink OdorValveOn
    0.75 0 0.75;...   %  Dark Pink OdorPokeOut
    0 0.75 0.75;...   %  Turquoise WaterPokeIn
    0 0 1;...         %  Navy Blue WaterValveOn
    0.75 0 0;...      %  Red WaterPokeOut
    0 1 0;...         %  Neon Green NextOdorPokeIn
    0.87 0.49 0];     %  Dark Purple Cue Signal
%% L & R Trials
% L trials
for stim_num = 1:num_stim
    
    % determine height of raster window
    num_trials = length(raster_events.L.trial_start{stim_num});
    raster_h = (num_trials / total_num_trials) * available_rasters_h;
    raster_u = raster_u - (raster_h + intrastim_space);
    
    trial_start_L = raster_events.L.trial_start{stim_num};
    align_event_time_L = raster_events.L.trial_events_time{stim_num}(align_ind, :);
    secondary_event_time_L = raster_events.L.trial_events_time{stim_num};
    
    raster_pos = [raster_l raster_u raster_w raster_h];
    
%     set(gcf, 'CurrentAxes', whole_fig); % re-reference whole figure
%     text(raster_l - textbuff, raster_u + raster_h/2, num2str(odorid(stim_num)));
    
    
    if num_trials > 0 % make sure there are trials of this type
        
%         axes('Position', raster_pos);
        
        [ref_spike_times.L{stim_num}, trial_inds.L{stim_num}, plot_handle] =...
            raster(spike_times, trial_start_L, align_event_time_L, window, secondary_event_time_L, 1);
        
%         set(gca, 'XTick', [], 'YTick', []);
        
        % Change colors of secondary event symbols
%         for sec_event_ind = 1:size(SECONDARY_EVENTS_COLORS, 1)
%             set(plot_handle.secondary_events(:, sec_event_ind), 'Color', SECONDARY_EVENTS_COLORS(sec_event_ind, :));
%         end
        
        % for first stim type, add the plot title
        %         if stim_num == 1
        %             f = findstr(spk_file, '\');
        %             title_str = [behav_file(end-14:end), ' - ', spk_file(f(end)+1:end)];
        %             title(ReplaceUnderscores(title_str));
        %         end
        
        % make a color bar showing choice
%         col = RASTER_COLORS.L;
%         choicebar_pos = [choicebar_l raster_u choicebar_w raster_h];
%         h = axes('Position', choicebar_pos);
%         set(h, 'Color', col, 'XColor', col, 'YColor', col);
%         set(h, 'YTick', [], 'XTick', []);
        
    else
        
        ref_spike_times.L{stim_num} = [];
        trial_inds.L{stim_num} = [];
        
    end
    
% R trials
    
    % determine height of raster window
    num_trials = length(raster_events.R.trial_start{stim_num});
    raster_h = (num_trials / total_num_trials) * available_rasters_h;
    raster_u = raster_u - (raster_h + intrachoice_space);
    
    trial_start_R = raster_events.R.trial_start{stim_num};
    align_event_time_R = raster_events.R.trial_events_time{stim_num}(align_ind, :);
    secondary_event_time_R = raster_events.R.trial_events_time{stim_num};
    
    raster_pos = [raster_l raster_u raster_w raster_h];
    
    if num_trials > 0 % make sure there are trials of this type
        
        % h = axes('Position', raster_pos);
%         axes('Position', raster_pos);
        
        [ref_spike_times.R{stim_num}, trial_inds.R{stim_num}, plot_handle] =...
            raster(spike_times, trial_start_R, align_event_time_R, window, secondary_event_time_R, 1);
        
%         set(gca, 'XTick', [], 'YTick', []);
        
        % Change colors of secondary event symbols
%         for sec_event_ind = 1:size(SECONDARY_EVENTS_COLORS, 1)
%             set(plot_handle.secondary_events(:, sec_event_ind), 'Color', SECONDARY_EVENTS_COLORS(sec_event_ind, :));
%         end
        
        % make a color bar showing choice
%         col = RASTER_COLORS.R;
%         choicebar_pos = [choicebar_l raster_u choicebar_w raster_h];
%         h = axes('Position', choicebar_pos);
%         set(h, 'Color', col, 'XColor', col, 'YColor', col);
%         set(h, 'YTick', [], 'XTick', []);
        
    else
        
        ref_spike_times.R{stim_num} = [];
        trial_inds.R{stim_num} = [];
        
    end
    
    
end

%% PSTHS

psth_l = raster_l;
psth_w = raster_w;
psth_h = 0.15;
psth_vertbuff = 0.01;
psth_u = raster_u - psth_h - psth_vertbuff;
psth_pos = [psth_l psth_u psth_w psth_h];


% axes('Position', psth_pos);
% 
% hold on;
% 
init_psth_window = window + [-((RESOLUTION/1000)/2) +((RESOLUTION/1000)/2)]; %b/c raster output has been padded

%% L & R Trials
for stim_num = 1:num_stim
    
    % L trials
    
%     if max(trial_inds.L{stim_num}) >= MIN_TRIAL_PLOTS % make sure there are enough trials of this type to plot
%         no_plot_flag = 0;
%     else
%         no_plot_flag = 1;
%     end
    
    [ref_psth.L.correct{stim_num}, plot_handle] = psth(ref_spike_times.L{stim_num}, trial_inds.L{stim_num},...
        init_psth_window, PSTH_SMOOTH_FACTOR, 1);
%     set(plot_handle, 'Color', RASTER_COLORS.L);
    
    
    % R trials
    
%     if max(trial_inds.R{stim_num}) >= MIN_TRIAL_PLOTS % make sure there are enough trials of this type to plot
%         no_plot_flag = 0;
%     else
%         no_plot_flag = 1;
%     end
    
    [ref_psth.R.correct{stim_num}, plot_handle] = psth(ref_spike_times.R{stim_num}, trial_inds.R{stim_num},...
        init_psth_window, PSTH_SMOOTH_FACTOR, 1);
%     set(plot_handle, 'Color', RASTER_COLORS.R);
    
    
%     xlabel(['Time relative to ', trial_events{align_ind}, ' (s)']);
%     ylabel('Firing rate (spks/s)');
    
end


%% package output structure

raster_info.ref_spike_times = ref_spike_times;
raster_info.trial_inds = trial_inds;
raster_info.raster_events = raster_events;
raster_info.stim_params = stim_params;
raster_info.align_ind = align_ind;
raster_info.window = window;
raster_info.beTI = before_trial_res;
raster_info.afTI = after_trial_res;



return