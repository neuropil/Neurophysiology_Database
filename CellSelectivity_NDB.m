%% CellSelectivity
% Calculate selectivity for a cell using ROC analysis, across multiple
% epochs. The selectivities currently calulated are:
% leftchoice_vs_rightchoice
% leftchoice_vs_rightchoice_5050trials
% leftchoice_vs_rightchoice_correct
% leftchoice_vs_rightchoice_error
% leftodors_vs_rightodors_pure - REMOVED
% leftchoice_vs_rightchoice_correct_pair1
% leftchoice_vs_rightchoice_correct_pair2
% leftchoice_vs_rightchoice_error_pair1
% leftchoice_vs_rightchoice_error_pair2
% odorA_vs_odorB_novel_correct
% odorA_vs_odorB_novel_error
% odorA_vs_odorB_novel_left
% odorA_vs_odorB_novel_right
% odorA50_vs_odorA100_left
% odorA50_vs_odorA100_right
% each odor vs every other odor
% pair1_vs_pair2_left_correct
% pair1_vs_pair2_right_correct
% rewarded_vs_unrewarded
% familiar_vs_novel_left_correct - REMOVED
% familiar_vs_novel_right_correct - REMOVED
%
% Added 3/19/09 and after:
% leftchoice_vs_rightchoice_easy
% leftchoice_vs_rightchoice_hard
% leftchoice_vs_rightchoice_short
% leftchoice_vs_rightchoice_long
% leftchoice_vs_rightchoice_pair1
% leftchoice_vs_rightchoice_pair2
%
% To add a selectivity calculation, add it to this file.
%
% USAGE: [selectivity, roc_p, activity_info] = CellSelectivity(raster_info, num_iter, analysis_to_perform)
% EXAMPLE: 
%
% INPUTS:  raster_info - structure created by batch_SC_raster_info
%              containing spike times relative to trial events, and
%              stimulus info
%          num_iter - # of iterations to calculate significance
%          analysis_to_perform (optional) - name of analysis to perform
%              (if not specified, all analyses are performed)
%          local_epochs (optional) - if specified, overrides the epochs
%              listed in GetSCGlobals
%
% OUTPUTS: selectivity - structure of selectivity values
%          roc_p - structure of significance (p) values for selectivity
%          activity_info - structure of additional activity information
%              .trial_spikes - spikes in each epoch of each trial
%              .stim_params - info about the stimulus from
%                  TrialClassification_Solo
%              .current_globals - the global variables used by this file
%                  (defined in GetSCGlobals.m)
%
% 9/17/07 - GF
%
% UPDATES:
% $ 12/27/07 - takes raster_info as input (spike times relative to trial
%            events already calculated and saved)
% $ 1/8/08 - only calculates selectivity in desired epoch, instead of all epochs
%          (saves a lot of time)
% $ 1/8/08 - Trials are excluded based on movement times that are too long
%          (since the rat may not have moved directly to the port)
% $ 2/7/08 - analysis_to_perform input added; allows user to specify only
%          1 analysis to be performed (reduces redundancy)
% $ 4/9/09 - added ability to define start AND stop of epochs w/ time
%              relative to trial events
% $ 7/2/11 - added option of inputting epochs, instead of using the epochs
%              in GetSCGlobals

function [selectivity, roc_p, activity_info] = CellSelectivity_NDB(raster_info, num_iter, analysis_to_perform, local_epochs)

if nargin < 3 % analysis_name not specified
    analysis_to_perform = 'all';
end

GetNeuroGlobals; % assign global variables
global MAX_MOVEMENT_TIME_TO_WATER; % max time from odorpokeout to waterpokein
% global MAX_MOVEMENT_TIME_TO_ODOR; % max time from waterpokeout to nextodorpokein
global epochs; % define epochs for selectivity analyses
global NOVEL_IND; %which mix_ind is the familiar and which is the novel odor
global GT_CUTOFF; % go tone delays shorter than this are "short" trials; longer are "long" trials

if nargin == 4 % use the epochs inputted, not the epochs in GetSCGlobals - 7/2/11
    epochs = local_epochs;
end

global EpochNames;
global EpochINFO;

% selectivity.EpochNames = EpochNames;
% selectivity.EpochDescription = EpochINFO;

%% Deconstruct the raster_info structure
num_stim = length(raster_info.raster_events.L.trial_start);
ref_spike_times = raster_info.ref_spike_times;
raster_events = raster_info.raster_events;
trial_inds = raster_info.trial_inds;
% stim_params = raster_info.stim_params;
align_ind = raster_info.align_ind;

%% For backward compatibility for reaction time cells, assign NaN to the gotone time for each trial
if size(raster_events.L.trial_events_time{1}, 1) == 7 % this is a rxn time cell
    for stim_num = 1:num_stim
        raster_events.L.trial_events_time{stim_num}(8, :) = NaN;
        raster_events.R.trial_events_time{stim_num}(8, :) = NaN;
    end
end

%%% Calculate the number of spikes in each epoch for each stimulus type %%%

for stim_num = 1:num_stim
    
    %% L trials

    if length(raster_events.L.trial_start{stim_num}) > 0 % make sure there are some trials of this type
        
        num_current_trials = length(raster_events.L.trial_start{stim_num});

        for trial_ind = 1:num_current_trials
            
            trial_spike_times = ref_spike_times.L{stim_num}(trial_inds.L{stim_num} == trial_ind);
            
            for epoch_ind = 1:length(epochs)
                
                %% get epoch start and stop times
                
                e = epochs{epoch_ind}; % current epoch

                % subtract the time of the event the rasters are aligned to (odorpokein)
                align_time = raster_events.L.trial_events_time{stim_num}(align_ind, trial_ind);
                    
                if length(e) == 2 % epoch boundaries correspond to trial events
                    
                    epoch_start = raster_events.L.trial_events_time{stim_num}(e(1), trial_ind) - align_time;
                    epoch_stop  = raster_events.L.trial_events_time{stim_num}(e(2), trial_ind) - align_time;

                elseif length(e) == 3 % one epoch boundaries does NOT correspond to a trial event
                                       
                    if e(1) == 888 % start boundary is determined by a time relative to a trial event
                        
                        epoch_stop = raster_events.L.trial_events_time{stim_num}(e(3), trial_ind) - align_time;
                        epoch_start = epoch_stop + e(2);
                        
                    elseif e(1) == 999 % stop boundary is determined by a time relative to a trial event
                        
                        epoch_start = raster_events.L.trial_events_time{stim_num}(e(2), trial_ind) - align_time;
                        epoch_stop = epoch_start + e(3);
                        
                    else
                        error('Flags for epoch start/stop other than trial events must be 888 or 999\n');
                    end
                    
                elseif length(e) == 6 % both epoch boundaries do NOT correspond to trial events
                        % format (e.g.): [888, time before event1, event1, 999, event2, time after event2] 
                    
                    % get start boundary    
                    if e(1) == 888 % start boundary is determined by a time BEFORE a trial event
                        
                        tmp = raster_events.L.trial_events_time{stim_num}(e(3), trial_ind) - align_time;
                        epoch_start = tmp + e(2);
                    
                    elseif e(1) == 999 % start boundary is determined by a time AFTER a trial event
                        
                        tmp = raster_events.L.trial_events_time{stim_num}(e(2), trial_ind) - align_time;
                        epoch_start = tmp + e(3);
                        
                    end
                        
                    % get stop boundary    
                    if e(4) == 888 % stop boundary is determined by a time BEFORE a trial event
                        
                        tmp = raster_events.L.trial_events_time{stim_num}(e(6), trial_ind) - align_time;
                        epoch_stop = tmp + e(5);
                    
                    elseif e(4) == 999 % stop boundary is determined by a time AFTER a trial event
                        
                        tmp = raster_events.L.trial_events_time{stim_num}(e(5), trial_ind) - align_time;
                        epoch_stop = tmp + e(6);
                        
                    end
                        
                else
                    error('each epoch entry must be a vector of length 2, 3, or 6.\n');
                end
                
                %% now that we have the epoch start and stop times, get the
                %% firing rate during the epoch
                
                trial_spikes.L{stim_num, epoch_ind}(trial_ind) = length(find((trial_spike_times > epoch_start) &...
                    (trial_spike_times <= epoch_stop))) / (epoch_stop - epoch_start);
                
            end
            
        end

        num_trials.L(stim_num) = num_current_trials;
                
    else % if no trials of this type
        
        for epoch_ind = 1:length(epochs)
            trial_spikes.L{stim_num, epoch_ind} = [];            
        end
        
        num_trials.L(stim_num) = 0;
    
    end

    %% R trials

    if length(raster_events.R.trial_start{stim_num}) > 0 % make sure there are some trials of this type
        
        num_current_trials = length(raster_events.R.trial_start{stim_num});

        for trial_ind = 1:num_current_trials
            
            trial_spike_times = ref_spike_times.R{stim_num}(trial_inds.R{stim_num} == trial_ind);
            
            for epoch_ind = 1:length(epochs)
                
                %% get epoch start and stop times
                
                e = epochs{epoch_ind}; % current epoch

                % subtract the time of the event the rasters are aligned to (odorpokein)
                align_time = raster_events.R.trial_events_time{stim_num}(align_ind, trial_ind);
                    
                if length(e) == 2 % epoch boundaries correspond to trial events
                    
                    epoch_start = raster_events.R.trial_events_time{stim_num}(e(1), trial_ind) - align_time;
                    epoch_stop  = raster_events.R.trial_events_time{stim_num}(e(2), trial_ind) - align_time;

                elseif length(e) == 3 % one or more epoch boundaries does NOT correspond to a trial event
                                       
                    if e(1) == 888 % start boundary is determined by a time relative to a trial event
                        
                        epoch_stop = raster_events.R.trial_events_time{stim_num}(e(3), trial_ind) - align_time;
                        epoch_start = epoch_stop + e(2);
                        
                    elseif e(1) == 999 % stop boundary is determined by a time relative to a trial event
                        
                        epoch_start = raster_events.R.trial_events_time{stim_num}(e(2), trial_ind) - align_time;
                        epoch_stop = epoch_start + e(3);
                        
                    else
                        error('Flags for epoch start/stop other than trial events must be 888 or 999\n');
                    end
                    
                elseif length(e) == 6 % both epoch boundaries do NOT correspond to trial events
                        % format (e.g.): [888, time before event1, event1, 999, event2, time after event2] 
                    
                    % get start boundary    
                    if e(1) == 888 % start boundary is determined by a time BEFORE a trial event
                        
                        tmp = raster_events.R.trial_events_time{stim_num}(e(3), trial_ind) - align_time;
                        epoch_start = tmp + e(2);
                    
                    elseif e(1) == 999 % start boundary is determined by a time AFTER a trial event
                        
                        tmp = raster_events.R.trial_events_time{stim_num}(e(2), trial_ind) - align_time;
                        epoch_start = tmp + e(3);
                        
                    end
                        
                    % get stop boundary    
                    if e(4) == 888 % stop boundary is determined by a time BEFORE a trial event
                        
                        tmp = raster_events.R.trial_events_time{stim_num}(e(6), trial_ind) - align_time;
                        epoch_stop = tmp + e(5);
                    
                    elseif e(4) == 999 % start boundary is determined by a time AFTER a trial event
                        
                        tmp = raster_events.R.trial_events_time{stim_num}(e(5), trial_ind) - align_time;
                        epoch_stop = tmp + e(6);
                        
                    end
                        
                else
                    error('each epoch entry must be a vector of length 2, 3, or 6.\n');
                end
                
                %% now that we have the epoch start and stop times, get the
                %% firing rate during the epoch
                
                trial_spikes.R{stim_num, epoch_ind}(trial_ind) = length(find((trial_spike_times > epoch_start) &...
                    (trial_spike_times <= epoch_stop))) / (epoch_stop - epoch_start);
                                
            end
            
        end

        num_trials.R(stim_num) = num_current_trials;
                
    else % if no trials of this type
        
        for epoch_ind = 1:length(epochs)
            trial_spikes.R{stim_num, epoch_ind} = [];            
        end
        
        num_trials.R(stim_num) = 0;
    
    end

end
    

%%% Calculate selectivity using ROC analysis %%%

%% Left choice vs. Right choice (all trial types): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

%     relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity
     relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity


    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = 1:num_stim

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = 1:num_stim

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            
            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% Left choice vs. Right choice (50/50 trials only): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_5050trials';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    %relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity
    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    stim_5050_inds = find(raster_info.stim_params.percent_A == 50); % find the 50/50 stimulus

    
    if length(stim_5050_inds) > 0 % if there are no 50/50 trials, skip this calculation

        for epoch_ind = relevant_epochs

            L_num_trials = []; % initialize
            R_num_trials = []; % initialize
            L_spikes = []; % initialize
            R_spikes = []; % initialize

            % L trials
            for stim_num = stim_5050_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.L.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                L_num_trials = [L_num_trials; length(valid_trials)];
                L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

            end

            % R trials
            for stim_num = stim_5050_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.R.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                R_num_trials = [R_num_trials; length(valid_trials)];
                R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

            end

            class_ID = [...
                -ones(sum(L_num_trials), 1);...
                 ones(sum(R_num_trials), 1);...
                ];

            observed_vals = [L_spikes; R_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else % set values to NaN if not applicable

        selectivity.(analysis_name) = NaN * ones(1, length(epochs));
        roc_p.(analysis_name) = NaN * ones(1, length(epochs));

    end

end


%% L vs. R - CORRECT trials only: -1=Left; 1=Right
%% 50-50 trials not included

analysis_name = 'leftchoice_vs_rightchoice_correct';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

%     relevant_epochs = [1:4]; % the epochs in which to calculate this selectivity
    relevant_epochs = 1:length(epochs);
    
    L_correct_stim_inds = find(raster_info.stim_params.reward_side == 1);
    R_correct_stim_inds = find(raster_info.stim_params.reward_side == 2);

    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = L_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = R_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% L vs. R - ERROR trials only: -1=Left; 1=Right
%% 50-50 trials not included

analysis_name = 'leftchoice_vs_rightchoice_error';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    L_error_stim_inds = find(raster_info.stim_params.reward_side == 2);
    R_error_stim_inds = find(raster_info.stim_params.reward_side == 1);

    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = L_error_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = R_error_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% L vs. R - BY ODOR PAIR, CORRECT trials only: -1=Left; 1=Right
%% 50-50 trials not included

for ind = 1:2 % 2 mix inds: familiar and novel, or mixtures and pure

    pair_num = ind;
    
    analysis_name = strcat('leftchoice_vs_rightchoice_correct_pair', num2str(pair_num));

    if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

        % initialize
        selectivity.(analysis_name) = NaN(1, length(epochs));
        roc_p.(analysis_name) = NaN(1, length(epochs));

        relevant_epochs = [2:4 9:11]; % the epochs in which to calculate this selectivity

        L_correct_stim_inds = find((raster_info.stim_params.reward_side == 1) & (raster_info.stim_params.mix_ind == pair_num));
        R_correct_stim_inds = find((raster_info.stim_params.reward_side == 2) & (raster_info.stim_params.mix_ind == pair_num));

        for epoch_ind = relevant_epochs

            L_num_trials = []; % initialize
            R_num_trials = []; % initialize
            L_spikes = []; % initialize
            R_spikes = []; % initialize

            % L trials
            for stim_num = L_correct_stim_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.L.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                L_num_trials = [L_num_trials; length(valid_trials)];
                L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

            end

            % R trials
            for stim_num = R_correct_stim_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.R.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                R_num_trials = [R_num_trials; length(valid_trials)];
                R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

            end

            class_ID = [...
                -ones(sum(L_num_trials), 1);...
                 ones(sum(R_num_trials), 1);...
                ];

            observed_vals = [L_spikes; R_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end

    end

end


%% L vs. R - BY ODOR PAIR, ERROR trials only: -1=Left; 1=Right
%% 50-50 trials not included

for ind = 1:2 % 2 mix inds: familiar and novel, or mixtures and pure

    pair_num = ind;
    
    analysis_name = strcat('leftchoice_vs_rightchoice_error_pair', num2str(pair_num));

    if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

        % initialize
        selectivity.(analysis_name) = NaN(1, length(epochs));
        roc_p.(analysis_name) = NaN(1, length(epochs));

        relevant_epochs = [2:4 9:11]; % the epochs in which to calculate this selectivity

        L_error_stim_inds = find((raster_info.stim_params.reward_side == 2) & (raster_info.stim_params.mix_ind == pair_num));
        R_error_stim_inds = find((raster_info.stim_params.reward_side == 1) & (raster_info.stim_params.mix_ind == pair_num));

        for epoch_ind = relevant_epochs

            L_num_trials = []; % initialize
            R_num_trials = []; % initialize
            L_spikes = []; % initialize
            R_spikes = []; % initialize

            % L trials
            for stim_num = L_error_stim_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.L.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                L_num_trials = [L_num_trials; length(valid_trials)];
                L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

            end

            % R trials
            for stim_num = R_error_stim_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.R.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                R_num_trials = [R_num_trials; length(valid_trials)];
                R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

            end

            class_ID = [...
                -ones(sum(L_num_trials), 1);...
                 ones(sum(R_num_trials), 1);...
                ];

            observed_vals = [L_spikes; R_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end

    end

end


%% Odor A vs. Odor B for the novel (and always pure) pair only - CORRECT trials only: -1=OdorA; 1=OdorB
%% Assume that correct for odor A is L, and correct for odor B is R

analysis_name = 'odorA_vs_odorB_novel_correct';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity

    if length(stim_params.mix_ind) == 4 % this is a 'novel odors' rat (i.e., no mixtures)

        odorA_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 100));
        odorB_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 0));

        for epoch_ind = relevant_epochs

            % Odor A trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{odorA_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            A_num_trials = length(valid_trials);
            A_spikes = trial_spikes.L{odorA_stim_num, epoch_ind}(valid_trials)';


            % Odor B trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{odorB_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            B_num_trials = length(valid_trials);
            B_spikes = trial_spikes.R{odorB_stim_num, epoch_ind}(valid_trials)';


            class_ID = [...
                -ones(sum(A_num_trials), 1);...
                 ones(sum(B_num_trials), 1);...
                ];

            observed_vals = [A_spikes; B_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else  % this is a 'mixtures' rat, so set selectivity and roc_p to NaN for all epochs

        selectivity.(analysis_name)(1:length(epochs)) = NaN;
        roc_p.(analysis_name)(1:length(epochs)) = NaN;

    end

end


%% Odor A vs. Odor B for the novel (and always pure) pair only - ERROR trials only: -1=OdorA; 1=OdorB
%% Assume that correct for odor A is L, and correct for odor B is R

analysis_name = 'odorA_vs_odorB_novel_error';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity

    if length(stim_params.mix_ind) == 4 % this is a 'novel odors' rat (i.e., no mixtures)

        odorA_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 100));
        odorB_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 0));

        for epoch_ind = relevant_epochs

            % Odor A trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{odorA_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            A_num_trials = length(valid_trials);
            A_spikes = trial_spikes.R{odorA_stim_num, epoch_ind}(valid_trials)';


            % Odor B trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{odorB_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            B_num_trials = length(valid_trials);
            B_spikes = trial_spikes.L{odorB_stim_num, epoch_ind}(valid_trials)';


            class_ID = [...
                -ones(sum(A_num_trials), 1);...
                 ones(sum(B_num_trials), 1);...
                ];

            observed_vals = [A_spikes; B_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else  % this is a 'mixtures' rat, so set selectivity and roc_p to NaN for all epochs

        selectivity.(analysis_name)(1:length(epochs)) = NaN;
        roc_p.(analysis_name)(1:length(epochs)) = NaN;

    end

end


%% Odor A vs. Odor B for the novel (and always pure) pair only - LEFT trials only: -1=OdorA; 1=OdorB
%% Assume that correct for odor A is L, and correct for odor B is R, thus...
%% This can also be thought of as novel-A-correct vs. novel-B-error.

analysis_name = 'odorA_vs_odorB_novel_left';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity

    if length(stim_params.mix_ind) == 4 % this is a 'novel odors' rat (i.e., no mixtures)

        odorA_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 100));
        odorB_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 0));

        for epoch_ind = relevant_epochs

            % Odor A trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{odorA_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            A_num_trials = length(valid_trials);
            A_spikes = trial_spikes.L{odorA_stim_num, epoch_ind}(valid_trials)';


            % Odor B trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{odorB_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            B_num_trials = length(valid_trials);
            B_spikes = trial_spikes.L{odorB_stim_num, epoch_ind}(valid_trials)';


            class_ID = [...
                -ones(sum(A_num_trials), 1);...
                 ones(sum(B_num_trials), 1);...
                ];

            observed_vals = [A_spikes; B_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else  % this is a 'mixtures' rat, so set selectivity and roc_p to NaN for all epochs

        selectivity.(analysis_name)(1:length(epochs)) = NaN;
        roc_p.(analysis_name)(1:length(epochs)) = NaN;

    end

end


%% Odor A vs. Odor B for the novel (and always pure) pair only - RIGHT trials only: -1=OdorA; 1=OdorB
%% Assume that correct for odor A is L, and correct for odor B is R, thus...
%% This can also be thought of as novel-A-error vs. novel-B-correct.

analysis_name = 'odorA_vs_odorB_novel_right';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity

    if length(stim_params.mix_ind) == 4 % this is a 'novel odors' rat (i.e., no mixtures)

        odorA_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 100));
        odorB_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 0));

        for epoch_ind = relevant_epochs

            % Odor A trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{odorA_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            A_num_trials = length(valid_trials);
            A_spikes = trial_spikes.R{odorA_stim_num, epoch_ind}(valid_trials)';


            % Odor B trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{odorB_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            B_num_trials = length(valid_trials);
            B_spikes = trial_spikes.R{odorB_stim_num, epoch_ind}(valid_trials)';


            class_ID = [...
                -ones(sum(A_num_trials), 1);...
                 ones(sum(B_num_trials), 1);...
                ];

            observed_vals = [A_spikes; B_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else  % this is a 'mixtures' rat, so set selectivity and roc_p to NaN for all epochs

        selectivity.(analysis_name)(1:length(epochs)) = NaN;
        roc_p.(analysis_name)(1:length(epochs)) = NaN;

    end

end


%% 50-50 mixtures vs. 100-0 mixtures - LEFT, 'correct' (for pure odor) trials only: -1=50-50; 1=100-0

analysis_name = 'odorA50_vs_odorA100_left';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity

    odorA50_stim_num = find(stim_params.percent_A == 50); % find the 50/50 stimulus
    if length(odorA50_stim_num) > 0  % find the 100/0 stimulus
        odorA100_stim_num = find((stim_params.percent_A == 100) & (stim_params.mix_ind == stim_params.mix_ind(odorA50_stim_num)));
    end

    if length(odorA50_stim_num) > 1

        error('There is more than one 50/50 stimulus type. Recode CellSelectivity.\n');

    elseif length(odorA50_stim_num) == 1 % if there are no 50/50 trials, skip this calculation

        for epoch_ind = relevant_epochs

            % 50/50 trials trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{odorA50_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            odorA50_num_trials = length(valid_trials);
            odorA50_spikes = trial_spikes.L{odorA50_stim_num, epoch_ind}(valid_trials)';

            % 100/0 trials trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{odorA100_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            odorA100_num_trials = length(valid_trials);
            odorA100_spikes = trial_spikes.L{odorA100_stim_num, epoch_ind}(valid_trials)';


            class_ID = [...
                -ones(sum(odorA50_num_trials), 1);...
                 ones(sum(odorA100_num_trials), 1);...
                ];

            observed_vals = [odorA50_spikes; odorA100_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else % set values to NaN if not applicable

        selectivity.(analysis_name) = NaN * ones(1, length(epochs));
        roc_p.(analysis_name) = NaN * ones(1, length(epochs));

    end

end


%% 50-50 mixtures vs. 100-0 mixtures - RIGHT, 'correct' (for pure odors) trials only: -1=50-50; 1=100-0

analysis_name = 'odorA50_vs_odorA100_right';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity

    odorA50_stim_num = find(stim_params.percent_A == 50); % find the 50/50 stimulus
    if length(odorA50_stim_num) > 0  % find the 100/0 stimulus
        odorA100_stim_num = find((stim_params.percent_A == 0) & (stim_params.mix_ind == stim_params.mix_ind(odorA50_stim_num)));
    end
    
    if length(odorA50_stim_num) > 1

        error('There is more than one 50/50 stimulus type. Recode CellSelectivity.\n');

    elseif length(odorA50_stim_num) == 1 % if there are no 50/50 trials, skip this calculation

        for epoch_ind = relevant_epochs

            % 50/50 trials trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{odorA50_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            odorA50_num_trials = length(valid_trials);
            odorA50_spikes = trial_spikes.R{odorA50_stim_num, epoch_ind}(valid_trials)';

            % 100/0 trials trials

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{odorA100_stim_num};
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            odorA100_num_trials = length(valid_trials);
            odorA100_spikes = trial_spikes.R{odorA100_stim_num, epoch_ind}(valid_trials)';


            class_ID = [...
                -ones(sum(odorA50_num_trials), 1);...
                 ones(sum(odorA100_num_trials), 1);...
                ];

            observed_vals = [odorA50_spikes; odorA100_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end   

    else % set values to NaN if not applicable

        selectivity.(analysis_name) = NaN * ones(1, length(epochs));
        roc_p.(analysis_name) = NaN * ones(1, length(epochs));

    end

end

%% Outcome: No reward vs. reward (all trial types): -1=no reward; 1=reward
% Note that for G017, trials can be correct w/o being rewarded!

analysis_name = 'rewarded_vs_unrewarded_no50';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = ([2,4,5,6,7]); % the epochs in which to calculate this selectivity
    %relevant_epochs = [1:length(epochs)]; % the epochs in which to calculate this selectivity
    left_odors = [2,4,6];
    right_odors = [1,3,5];

    % count the number of rewarded and unrewarded trials (based on whether the
    % watervalve went on)
    num_rewarded = []; % initialize
    num_unrewarded = []; % initialize
    for trial_type_ind = 1:length(left_odors) % L trials
        num_rewarded = [num_rewarded length(raster_info.raster_events.L.trial_start{left_odors(trial_type_ind)})];
        num_unrewarded = [num_unrewarded length(raster_info.raster_events.L.trial_start{right_odors(trial_type_ind)})];
    end
    for trial_type_ind = 1:length(right_odors) % R trials
        num_rewarded = [num_rewarded length(raster_info.raster_events.R.trial_start{right_odors(trial_type_ind)})];
        num_unrewarded = [num_unrewarded length(raster_info.raster_events.R.trial_start{left_odors(trial_type_ind)})];
    end

    class_ID = [...
        -ones(sum(num_unrewarded), 1);...
         ones(sum(num_rewarded), 1);...
        ];


    for epoch_ind = relevant_epochs

        rewarded_spikes = []; % initialize
        unrewarded_spikes = []; % initialize

        % L trials
        for stim_num = 1:length(left_odors)

            rewarded_spikes = [rewarded_spikes; trial_spikes.L{left_odors(stim_num), epoch_ind}'];
            unrewarded_spikes = [unrewarded_spikes; trial_spikes.L{right_odors(stim_num), epoch_ind}'];

        end

        % R trials
        for stim_num = 1:length(right_odors)

            rewarded_spikes = [rewarded_spikes; trial_spikes.R{right_odors(stim_num), epoch_ind}'];
            unrewarded_spikes = [unrewarded_spikes; trial_spikes.R{left_odors(stim_num), epoch_ind}'];

        end

        observed_vals = [unrewarded_spikes; rewarded_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;
        
        %% save spikerates in rewarded and unrewarded trials
        outcome_info(epoch_ind).rewarded_spikerate = rewarded_spikes;
        outcome_info(epoch_ind).unrewarded_spikerate = unrewarded_spikes;

    end

end

% %% One odor vs. all other odors: -1=other; 1=odor of interest
% 
% relevant_epochs = [2:3]; % the epochs in which to calculate this selectivity
% 
% pure_odors = find((stim_params.percent_A == 0) | (stim_params.percent_A == 100));
% num_odors = length(pure_odors);
% 
% for odor_ind = 1:num_odors
%     
%     analysis_name = strcat('odor', num2str(odor_ind));
%     
%     if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified
%     
%         % initialize
%         selectivity.(analysis_name) = NaN(1, length(epochs));
%         roc_p.(analysis_name) = NaN(1, length(epochs));
% 
%         odor_of_interest_stim_num = pure_odors(odor_ind);
%         other_odors_stim_num = pure_odors(pure_odors ~= odor_of_interest_stim_num);
% 
%         class_ID = [...
%             -ones(sum(num_trials.L(other_odors_stim_num)), 1);...
%             -ones(sum(num_trials.R(other_odors_stim_num)), 1);...
%              ones(sum(num_trials.L(odor_of_interest_stim_num)), 1);...
%              ones(sum(num_trials.R(odor_of_interest_stim_num)), 1);...
%             ];
% 
%         for epoch_ind = relevant_epochs
% 
%             observed_vals = [...
%                 cell2mat(trial_spikes.L(other_odors_stim_num, epoch_ind)')';...
%                 cell2mat(trial_spikes.R(other_odors_stim_num, epoch_ind)')';...
%                 cell2mat(trial_spikes.L(odor_of_interest_stim_num, epoch_ind)')';...
%                 cell2mat(trial_spikes.R(odor_of_interest_stim_num, epoch_ind)')';...
%                 ];
% 
%             [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction
% 
%             selectivity.(analysis_name)(epoch_ind) = sel;
%             roc_p.(analysis_name)(epoch_ind) = rp;
% 
%         end
% 
%     end
%     
% end


%% Pair 1 L odor vs. Pair 2 L odor - CORRECT trials only: -1=pair1; 1=pair2
%% Assume that correct for odor A is L, and correct for odor B is R
%% NOTE!: This includes the pure odors for mixtures rats AND the
%% novel/familiar odors for novel odors rats.
%% Pair 1 represents the familiar pair for novel odors rats, and the mixtures pair for mixtures rats. 

analysis_name = 'pair1_vs_pair2_left_correct';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:4 9]; % the epochs in which to calculate this selectivity

    pair1_stim_num = find((stim_params.mix_ind == 1) & (stim_params.percent_A == 100));
    pair2_stim_num = find((stim_params.mix_ind == 2) & (stim_params.percent_A == 100));

    if length(pair2_stim_num) > 0 % only do this analysis if there are 2 pairs of odors
        
        for epoch_ind = relevant_epochs

            pair1_num_trials = []; % initialize
            pair2_num_trials = []; % initialize
            pair1_spikes = []; % initialize
            pair2_spikes = []; % initialize

            % pair1 trials
            for stim_ind = 1:length(pair1_stim_num)

                % only include certain trials, based on movement time
                r = raster_info.raster_events.L.trial_events_time{pair1_stim_num(stim_ind)};
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                pair1_num_trials = [pair1_num_trials; length(valid_trials)];
                pair1_spikes = [pair1_spikes; trial_spikes.L{pair1_stim_num(stim_ind), epoch_ind}(valid_trials)'];

            end
            
            % pair2 trials
            for stim_ind = 1:length(pair2_stim_num)

                % only include certain trials, based on movement time
                r = raster_info.raster_events.L.trial_events_time{pair2_stim_num(stim_ind)};
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                pair2_num_trials = [pair2_num_trials; length(valid_trials)];
                pair2_spikes = [pair2_spikes; trial_spikes.L{pair2_stim_num(stim_ind), epoch_ind}(valid_trials)'];

            end
            
            class_ID = [...
                -ones(sum(pair1_num_trials), 1);...
                 ones(sum(pair2_num_trials), 1);...
                ];

            observed_vals = [pair1_spikes; pair2_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end

    else  % only 1 odor pair was used in this session, so set selectivity and roc_p to NaN for all epochs

        selectivity.(analysis_name)(1:length(epochs)) = NaN;
        roc_p.(analysis_name)(1:length(epochs)) = NaN;

    end

end


%% Pair 1 R odor vs. Pair 2 R odor - CORRECT trials only: -1=pair1; 1=pair2
%% Assume that correct for odor A is L, and correct for odor B is R
%% NOTE!: This includes the pure odors for mixtures rats AND the
%% novel/familiar odors for novel odors rats.
%% Pair 1 represents the familiar pair for novel odors rats, and the mixtures pair for mixtures rats. 

analysis_name = 'pair1_vs_pair2_right_correct';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = [2:4 9]; % the epochs in which to calculate this selectivity

    pair1_stim_num = find((stim_params.mix_ind == 1) & (stim_params.percent_A == 0));
    pair2_stim_num = find((stim_params.mix_ind == 2) & (stim_params.percent_A == 0));

    if length(pair2_stim_num) > 0 % only do this analysis if there are 2 pairs of odors
        
        for epoch_ind = relevant_epochs

            pair1_num_trials = []; % initialize
            pair2_num_trials = []; % initialize
            pair1_spikes = []; % initialize
            pair2_spikes = []; % initialize

            % pair1 trials
            for stim_ind = 1:length(pair1_stim_num)

                % only include certain trials, based on movement time
                r = raster_info.raster_events.R.trial_events_time{pair1_stim_num(stim_ind)};
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                pair1_num_trials = [pair1_num_trials; length(valid_trials)];
                pair1_spikes = [pair1_spikes; trial_spikes.R{pair1_stim_num(stim_ind), epoch_ind}(valid_trials)'];

            end
            
            % pair2 trials
            for stim_ind = 1:length(pair2_stim_num)

                % only include certain trials, based on movement time
                r = raster_info.raster_events.R.trial_events_time{pair2_stim_num(stim_ind)};
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                pair2_num_trials = [pair2_num_trials; length(valid_trials)];
                pair2_spikes = [pair2_spikes; trial_spikes.R{pair2_stim_num(stim_ind), epoch_ind}(valid_trials)'];

            end
            
            class_ID = [...
                -ones(sum(pair1_num_trials), 1);...
                 ones(sum(pair2_num_trials), 1);...
                ];

            observed_vals = [pair1_spikes; pair2_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end

    else  % only 1 odor pair was used in this session, so set selectivity and roc_p to NaN for all epochs

        selectivity.(analysis_name)(1:length(epochs)) = NaN;
        roc_p.(analysis_name)(1:length(epochs)) = NaN;

    end

end


%% Left choice vs. Right choice: easy trials only (>=80% of either odor - i.e., [100/0, 80/20, 20/80, 0/100]): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_easy';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    stim_to_use = find_any(raster_info.stim_params.percent_A, [95 -95 80 -80], 'logical'); % find the 'easy' stimuli
    stim_to_use = find(stim_to_use == 1);
   
    
    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = stim_to_use

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = stim_to_use

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% Left choice vs. Right choice: hard trials only (<=60% of either odor - i.e., [60/40, 50/50, 40/60]): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_hard';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

% initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    stim_to_use = find_any(raster_info.stim_params.percent_A, [60 -60 50], 'logical'); % find the 'hard' stimuli
    stim_to_use = find(stim_to_use == 1);

    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = stim_to_use

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = stim_to_use

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% Left choice vs. Right short: short delay trials only

analysis_name = 'leftchoice_vs_rightchoice_short';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for L_trial_stims = 1:length(raster_info.raster_events.L.long_index)

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{L_trial_stims};
            s = raster_info.raster_events.L.long_index{L_trial_stims};
            
            if isempty(s)
                valid_trials = [];
            else
                valid_trials = find((r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER) & ~s);
            end
                    
            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{L_trial_stims, epoch_ind}(valid_trials)'];
            
        end

        % R trials
        for R_trial_stims = 1:length(raster_info.raster_events.R.long_index)

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{R_trial_stims};
            s = raster_info.raster_events.R.long_index{R_trial_stims};
            
            if isempty(s)
                valid_trials = [];
            else
                valid_trials = find((r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER) & ~s);
            end

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{R_trial_stims, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% Left choice vs. Right short: long delay trials only

analysis_name = 'leftchoice_vs_rightchoice_long';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for L_trial_stims = 1:length(raster_info.raster_events.L.long_index)

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{L_trial_stims};
            l = raster_info.raster_events.L.long_index{L_trial_stims};
            
            if isempty(l)
                valid_trials = [];
            else
                valid_trials = find((r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER) & l);
            end
                    
            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{L_trial_stims, epoch_ind}(valid_trials)'];
            
        end

        % R trials
        for R_trial_stims = 1:length(raster_info.raster_events.R.long_index)

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{R_trial_stims};
            l = raster_info.raster_events.R.long_index{R_trial_stims};
            
            if isempty(l)
                valid_trials = [];
            else
                valid_trials = find((r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER) & l);
            end

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{R_trial_stims, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end


%% L vs. R - BY ODOR PAIR, all trials (correct and error): -1=Left; 1=Right

for ind = 1:2 % 2 mix inds: familiar and novel, or mixtures and pure

    pair_num = ind;
    
    analysis_name = strcat('leftchoice_vs_rightchoice_pair', num2str(pair_num));

    if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

        % initialize
        selectivity.(analysis_name) = NaN(1, length(epochs));
        roc_p.(analysis_name) = NaN(1, length(epochs));

        relevant_epochs = [2:4]; % the epochs in which to calculate this selectivity

        stim_inds = find(raster_info.stim_params.mix_ind == pair_num);

        for epoch_ind = relevant_epochs

            L_num_trials = []; % initialize
            R_num_trials = []; % initialize
            L_spikes = []; % initialize
            R_spikes = []; % initialize

            % L trials
            for stim_num = stim_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.L.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                L_num_trials = [L_num_trials; length(valid_trials)];
                L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

            end

            % R trials
            for stim_num = stim_inds

                % only include certain trials, based on movement time
                r = raster_info.raster_events.R.trial_events_time{stim_num};

%                 if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                     valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%                 else % epochs related to movement towards water port
%                     valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%                 end
                
                valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

                R_num_trials = [R_num_trials; length(valid_trials)];
                R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

            end

            class_ID = [...
                -ones(sum(L_num_trials), 1);...
                 ones(sum(R_num_trials), 1);...
                ];

            observed_vals = [L_spikes; R_spikes];

            [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

            selectivity.(analysis_name)(epoch_ind) = sel;
            roc_p.(analysis_name)(epoch_ind) = rp;

        end

    end

end



%% Left choice vs. Right choice: CORRECT, easy trials only (>=80% of either odor - i.e., [100/0, 80/20, 20/80, 0/100]): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_correct_easy';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    L_correct_stim_inds = find_any(raster_info.stim_params.percent_A, [-95 -80], 'logical') &...
        ((raster_info.stim_params.reward_side == 1) | (raster_info.stim_params.reward_side == 1.5)); % find the correct, L, 'easy' stimuli
    R_correct_stim_inds = find_any(raster_info.stim_params.percent_A, [95 80], 'logical') &...
        ((raster_info.stim_params.reward_side == 2) | (raster_info.stim_params.reward_side == 1.5)); % find the correct, R, 'easy' stimuli

    L_correct_stim_inds = find(L_correct_stim_inds == 1);
    R_correct_stim_inds = find(R_correct_stim_inds == 1);
    
    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = L_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = R_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end

%% Left choice vs. Right choice: ERROR, easy trials only (>=80% of either odor - i.e., [100/0, 80/20, 20/80, 0/100]): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_error_easy';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    L_error_stim_inds = find_any(raster_info.stim_params.percent_A, [95 80], 'logical');
    R_error_stim_inds = find_any(raster_info.stim_params.percent_A, [-95 -80], 'logical');

    L_error_stim_inds = find(L_error_stim_inds == 1);
    R_error_stim_inds = find(R_error_stim_inds == 1);
    
    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = L_error_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = R_error_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end
%% Left choice vs. Right choice: CORRECT, hard trials only (<=60% of either odor - i.e., [60/40, 50/50, 40/60]): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_correct_hard';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    L_correct_stim_inds = find_any(raster_info.stim_params.percent_A, [-60 50], 'logical') &...
        ((raster_info.stim_params.reward_side == 1) | (raster_info.stim_params.reward_side == 1.5)); % find the correct, L, 'hard' stimuli
    R_correct_stim_inds = find_any(raster_info.stim_params.percent_A, [60 50], 'logical') &...
        ((raster_info.stim_params.reward_side == 2) | (raster_info.stim_params.reward_side == 1.5)); % find the correct, R, 'hard' stimuli

    L_correct_stim_inds = find(L_correct_stim_inds == 1);
    R_correct_stim_inds = find(R_correct_stim_inds == 1);
    
    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = L_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = R_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end

%% Left choice vs. Right choice: CORRECT, hard trials only (<=60% of either odor - i.e., [60/40, 50/50, 40/60]): -1=L; 1=R

analysis_name = 'leftchoice_vs_rightchoice_error_hard';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = 1:length(epochs); % the epochs in which to calculate this selectivity

    L_correct_stim_inds = find_any(raster_info.stim_params.percent_A, [60 50], 'logical');
    R_correct_stim_inds = find_any(raster_info.stim_params.percent_A, [-60 50], 'logical');

    L_correct_stim_inds = find(L_correct_stim_inds == 1);
    R_correct_stim_inds = find(R_correct_stim_inds == 1);
    
    for epoch_ind = relevant_epochs

        L_num_trials = []; % initialize
        R_num_trials = []; % initialize
        L_spikes = []; % initialize
        R_spikes = []; % initialize

        % L trials
        for stim_num = L_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.L.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            L_num_trials = [L_num_trials; length(valid_trials)];
            L_spikes = [L_spikes; trial_spikes.L{stim_num, epoch_ind}(valid_trials)'];

        end

        % R trials
        for stim_num = R_correct_stim_inds

            % only include certain trials, based on movement time
            r = raster_info.raster_events.R.trial_events_time{stim_num};

%             if (epoch_ind == 7 | epoch_ind == 8) % epochs related to reinitiation movement
%                 valid_trials = find(r(7, :) - r(6, :) <= MAX_MOVEMENT_TIME_TO_ODOR);
%             else % epochs related to movement towards water port
%                 valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
%             end
            
            valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);

            R_num_trials = [R_num_trials; length(valid_trials)];
            R_spikes = [R_spikes; trial_spikes.R{stim_num, epoch_ind}(valid_trials)'];

        end

        class_ID = [...
            -ones(sum(L_num_trials), 1);...
             ones(sum(R_num_trials), 1);...
            ];

        observed_vals = [L_spikes; R_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;

    end

end



%% Outcome: No reward vs. reward (all trial types): -1=no reward; 1=reward
% Note that for G017, trials can be correct w/o being rewarded!

analysis_name = 'rewarded_vs_unrewarded';

if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified

    % initialize
    selectivity.(analysis_name) = NaN(1, length(epochs));
    roc_p.(analysis_name) = NaN(1, length(epochs));

    relevant_epochs = (5:6); % the epochs in which to calculate this selectivity
    %relevant_epochs = [1:length(epochs)]; % the epochs in which to calculate this selectivity

    % count the number of rewarded and unrewarded trials (based on whether the
    % watervalve went on)
    num_rewarded = []; % initialize
    num_unrewarded = []; % initialize
    for trial_type_ind = 1:length(raster_info.raster_events.L.trial_start) % L trials
        num_rewarded = [num_rewarded sum(~isnan(raster_info.raster_events.L.trial_events_time{trial_type_ind}(5, :)))];
        num_unrewarded = [num_unrewarded sum(isnan(raster_info.raster_events.L.trial_events_time{trial_type_ind}(5, :)))];
    end
    for trial_type_ind = 1:length(raster_info.raster_events.R.trial_start) % R trials
        num_rewarded = [num_rewarded sum(~isnan(raster_info.raster_events.R.trial_events_time{trial_type_ind}(5, :)))];
        num_unrewarded = [num_unrewarded sum(isnan(raster_info.raster_events.R.trial_events_time{trial_type_ind}(5, :)))];
    end

    class_ID = [...
        -ones(sum(num_unrewarded), 1);...
         ones(sum(num_rewarded), 1);...
        ];


    for epoch_ind = relevant_epochs

        rewarded_spikes = []; % initialize
        unrewarded_spikes = []; % initialize

        % L trials
        for stim_num = 1:num_stim

            rewarded_trial_nums = find(~isnan(raster_info.raster_events.L.trial_events_time{stim_num}(5, :)));
            unrewarded_trial_nums = find(isnan(raster_info.raster_events.L.trial_events_time{stim_num}(5, :)));

            rewarded_spikes = [rewarded_spikes; trial_spikes.L{stim_num, epoch_ind}(rewarded_trial_nums)'];
            unrewarded_spikes = [unrewarded_spikes; trial_spikes.L{stim_num, epoch_ind}(unrewarded_trial_nums)'];

        end

        % R trials
        for stim_num = 1:num_stim

            rewarded_trial_nums = find(~isnan(raster_info.raster_events.R.trial_events_time{stim_num}(5, :)));
            unrewarded_trial_nums = find(isnan(raster_info.raster_events.R.trial_events_time{stim_num}(5, :)));

            rewarded_spikes = [rewarded_spikes; trial_spikes.R{stim_num, epoch_ind}(rewarded_trial_nums)'];
            unrewarded_spikes = [unrewarded_spikes; trial_spikes.R{stim_num, epoch_ind}(unrewarded_trial_nums)'];

        end

        observed_vals = [unrewarded_spikes; rewarded_spikes];

        [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction

        selectivity.(analysis_name)(epoch_ind) = sel;
        roc_p.(analysis_name)(epoch_ind) = rp;
        
        %% save spikerates in rewarded and unrewarded trials
        outcome_info(epoch_ind).rewarded_spikerate = rewarded_spikes;
        outcome_info(epoch_ind).unrewarded_spikerate = unrewarded_spikes;

    end

end


%% save additional information in activity_info structure
activity_info.trial_spikes = trial_spikes;
%activity_info.stim_params = stim_params;

% spike info specific to outcome analyses should be saved as well (so that
% we have access to spikerates when compiling analyses)
if exist('outcome_info');
    activity_info.outcome_info = outcome_info;
end

% save the global variables that this file used
global_struct = whos('global');
for i = 1:length(global_struct)
    if exist(global_struct(i).name)
        current_globals.(global_struct(i).name) = eval(global_struct(i).name);
    end
end

activity_info.current_globals = current_globals;

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalculateROC subfunction - gets the selectivity and significance,
%% given the observed_vals and class_ID (i.e., does the repetitive
%% heavy lifting for all the different selectivity measures)

function [selectivity, significance] = CalculateROC(observed_vals, class_ID, num_iter)

if length(observed_vals) ~= length(class_ID)
    error('observed_vals and class_ID must have the same length!');
end

if sum(isnan(observed_vals)) == length(observed_vals) % all observed_vals are NaN, so don't calculate ROC
    selectivity = NaN;
    significance = NaN;
    return;
end

ALPHA = 0.05; % not actually used

roc_area = roc_sv(observed_vals, class_ID, 'nofigure');

selectivity = 2 * (roc_area - 0.5);

% calculate significance of roc area

% if roc_area < 0.5, flip class_IDs
if roc_area < 0.5
    class_ID_for_sig = -class_ID;
else
    class_ID_for_sig = class_ID;            
end

if ((num_iter > 0) && (sum(class_ID == -1) > 1) && (sum(class_ID == 1) > 1)) % make sure there are >1 of each observation
    
    % bootstrap ROC -- from RE Strauss library
    [conf_int, significance] = RES_Bootstrp('roc_bootstrp',[0 1 0 1], num_iter, ALPHA, observed_vals, class_ID_for_sig, 0);

else
    
    significance = NaN;
    
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unperformed analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Left odors vs. Right odors (pure odor trials only): -1=L; 1=R
% 
% analysis_name = 'leftodors_vs_rightodors_pure';
% 
% if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified
% 
%     % initialize
%     selectivity.(analysis_name) = NaN(1, length(epochs));
%     roc_p.(analysis_name) = NaN(1, length(epochs));
% 
%     relevant_epochs = [2:3]; % the epochs in which to calculate this selectivity
% 
%     left_odors_stim_num = find((stim_params.percent_A == 100) & (stim_params.reward_side == 1));
%     right_odors_stim_num = find((stim_params.percent_A == 0) & (stim_params.reward_side == 2));
% 
%     class_ID = [...
%         -ones(sum(num_trials.L(left_odors_stim_num)), 1);...
%         -ones(sum(num_trials.R(left_odors_stim_num)), 1);...
%          ones(sum(num_trials.L(right_odors_stim_num)), 1);...
%          ones(sum(num_trials.R(right_odors_stim_num)), 1);...
%         ];
% 
%     for epoch_ind = relevant_epochs
% 
%         observed_vals = [...
%             cell2mat(trial_spikes.L(left_odors_stim_num, epoch_ind)')';...
%             cell2mat(trial_spikes.R(left_odors_stim_num, epoch_ind)')';...
%             cell2mat(trial_spikes.L(right_odors_stim_num, epoch_ind)')';...
%             cell2mat(trial_spikes.R(right_odors_stim_num, epoch_ind)')';...
%             ];
% 
%         [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction
% 
%         selectivity.(analysis_name)(epoch_ind) = sel;
%         roc_p.(analysis_name)(epoch_ind) = rp;
% 
%     end   
% 
% end
% 
% 

% 
% 
% %% familiar vs. novel left odor - CORRECT trials only: -1=familiar; 1=novel
% %% Assume that correct for odor A is L, and correct for odor B is R
% 
% analysis_name = 'familiar_vs_novel_left_correct';
% 
% if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified
% 
%     % initialize
%     selectivity.(analysis_name) = NaN(1, length(epochs));
%     roc_p.(analysis_name) = NaN(1, length(epochs));
% 
%     relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity
% 
%     if length(stim_params.mix_ind) == 4 % this is a 'novel odors' rat (i.e., no mixtures)
% 
%         familiar_stim_num = find((stim_params.mix_ind == FAMILIAR_IND) & (stim_params.percent_A == 100));
%         novel_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 100));
% 
%         for epoch_ind = relevant_epochs
% 
%             % familiar trials
% 
%             % only include certain trials, based on movement time
%             r = raster_info.raster_events.L.trial_events_time{familiar_stim_num};
%             valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
% 
%             familiar_num_trials = length(valid_trials);
%             familiar_spikes = trial_spikes.L{familiar_stim_num, epoch_ind}(valid_trials)';
% 
% 
%             % novel trials
% 
%             % only include certain trials, based on movement time
%             r = raster_info.raster_events.L.trial_events_time{novel_stim_num};
%             valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
% 
%             novel_num_trials = length(valid_trials);
%             novel_spikes = trial_spikes.L{novel_stim_num, epoch_ind}(valid_trials)';
% 
% 
%             class_ID = [...
%                 -ones(sum(familiar_num_trials), 1);...
%                  ones(sum(novel_num_trials), 1);...
%                 ];
% 
%             observed_vals = [familiar_spikes; novel_spikes];
% 
%             [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction
% 
%             selectivity.(analysis_name)(epoch_ind) = sel;
%             roc_p.(analysis_name)(epoch_ind) = rp;
% 
%         end   
% 
%     else  % this is a 'mixtures' rat, so set selectivity and roc_p to NaN for all epochs
% 
%         selectivity.(analysis_name)(1:length(epochs)) = NaN;
%         roc_p.(analysis_name)(1:length(epochs)) = NaN;
% 
%     end
% 
% end
% 
% 
% %% familiar vs. novel right odor - CORRECT trials only: -1=familiar; 1=novel
% %% Assume that correct for odor A is L, and correct for odor B is R
% 
% analysis_name = 'familiar_vs_novel_right_correct';
% 
% if (strcmpi(analysis_to_perform, analysis_name) || strcmpi(analysis_to_perform, 'all')) % only perform this analysis if specified
% 
%     % initialize
%     selectivity.(analysis_name) = NaN(1, length(epochs));
%     roc_p.(analysis_name) = NaN(1, length(epochs));
% 
%     relevant_epochs = [2:6]; % the epochs in which to calculate this selectivity
% 
%     if length(stim_params.mix_ind) == 4 % this is a 'novel odors' rat (i.e., no mixtures)
% 
%         familiar_stim_num = find((stim_params.mix_ind == FAMILIAR_IND) & (stim_params.percent_A == 0));
%         novel_stim_num = find((stim_params.mix_ind == NOVEL_IND) & (stim_params.percent_A == 0));
% 
%         for epoch_ind = relevant_epochs
% 
%             % familiar trials
% 
%             % only include certain trials, based on movement time
%             r = raster_info.raster_events.R.trial_events_time{familiar_stim_num};
%             valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
% 
%             familiar_num_trials = length(valid_trials);
%             familiar_spikes = trial_spikes.R{familiar_stim_num, epoch_ind}(valid_trials)';
% 
% 
%             % novel trials
% 
%             % only include certain trials, based on movement time
%             r = raster_info.raster_events.R.trial_events_time{novel_stim_num};
%             valid_trials = find(r(4, :) - r(3, :) <= MAX_MOVEMENT_TIME_TO_WATER);
% 
%             novel_num_trials = length(valid_trials);
%             novel_spikes = trial_spikes.R{novel_stim_num, epoch_ind}(valid_trials)';
% 
% 
%             class_ID = [...
%                 -ones(sum(familiar_num_trials), 1);...
%                  ones(sum(novel_num_trials), 1);...
%                 ];
% 
%             observed_vals = [familiar_spikes; novel_spikes];
% 
%             [sel, rp] = CalculateROC(observed_vals, class_ID, num_iter); % Calls subfunction
% 
%             selectivity.(analysis_name)(epoch_ind) = sel;
%             roc_p.(analysis_name)(epoch_ind) = rp;
% 
%         end   
% 
%     else  % this is a 'mixtures' rat, so set selectivity and roc_p to NaN for all epochs
% 
%         selectivity.(analysis_name)(1:length(epochs)) = NaN;
%         roc_p.(analysis_name)(1:length(epochs)) = NaN;
% 
%     end
% 
% end



