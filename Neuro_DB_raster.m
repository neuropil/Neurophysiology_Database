%% raster.m
% Makes a raster plot given spike and referencing event times. Plotting
% window should be set up externally.
%
% USAGE: [ref_spike_times, trial_inds, plot_handle] = raster(spike_times,
%         trial_start_times, event_times, window, secondary_event_times, 
%         no_plot_flag)
% EXAMPLE: [spk_times, trials_of_spikes, h_raster] = raster(spk_times, 
%           TrialStart, WaterValveOn, [-2 2]);
%
% INPUTS:  spike_times - times of spikes (IN SECONDS)
%          trial_start_times - 1xN vector of trial start stimes (IN SECONDS) (often
%              saved as TrialStart from behavioral protocols); N = number
%              of trials
%          event_times - 1xN vector of times of events to reference spikes to (e.g.,
%              WaterPokeIn), relative to trial start (IN SECONDS)
%
% OUTPUTS: ref_spike_times - spike times referenced to trial event
%          trial_inds - trials in which corresponding ref_spike_times
%          occured


function [ref_spike_times, trial_inds, ref_times] = Neuro_DB_raster(spike_times, trial_start_times, event_times)

GetSCGlobals;
global RESOLUTION;
%RESOLUTION = 1000; % plot spike times to 1 msec resolution - DO NOT CHANGE

window = [-3 3];

%% expand the window on either side to account for smoothing performed by psth
% plotting
window(1) = window(1) - ((RESOLUTION/1000)/2);
window(2) = window(2) + ((RESOLUTION/1000)/2);

% initialize
ref_spike_times = [];
trial_inds = [];
ref_times = [];

for trial_num = 1:length(trial_start_times)
    
    ref_time = trial_start_times(trial_num) + event_times(trial_num); % absolute time of event to reference spikes to
    
    if ~isnan(ref_time) % ref_time will be NaN if the particular event did not occur in this trial
        
        % trial_spike_times is a temporary variable
        trial_spike_times = spike_times(find((spike_times >= (ref_time + window(1))) &...
            (spike_times <= (ref_time + window(2))))) - ref_time;
        
        ref_spike_times = [ref_spike_times; trial_spike_times(:)];
        trial_inds = [trial_inds; (ones(length(trial_spike_times), 1) * trial_num)];
        ref_times = [ref_times; (ones(length(trial_spike_times), 1) * ref_time)];
        
    end
    
end

return