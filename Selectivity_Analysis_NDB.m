function [selectivity_info] = Selectivity_Analysis_NDB(raster_info)

ANALYSIS = {'leftchoice_vs_rightchoice','leftchoice_vs_rightchoice_correct',...
    'leftchoice_vs_rightchoice_error', 'leftchoice_vs_rightchoice_correct_easy',...
    'leftchoice_vs_rightchoice_correct_hard', 'leftchoice_vs_rightchoice_error_easy',...
    'leftchoice_vs_rightchoice_error_hard', 'rewarded_vs_unrewarded',...
    'leftchoice_vs_rightchoice_5050trials', 'rewarded_vs_unrewarded_no50',...
    'leftchoice_vs_rightchoice_short','leftchoice_vs_rightchoice_long'};

NUM_ITER = 500; % how many iterations to use for the bootstrap estimation of ROC significance 

for analysis = 1:length(ANALYSIS);
                
        % Call a separate function to calculate selectivity for the cell
        [selectivity, roc_p, activity_info] = CellSelectivity_NDB(raster_info, NUM_ITER, char(ANALYSIS{analysis}));
        
        % package structure to save
        selectivity_info.selectivity.(char(ANALYSIS{analysis})) = selectivity.(char(ANALYSIS{analysis}));
        selectivity_info.roc_p.(char(ANALYSIS{analysis})) =  roc_p.(char(ANALYSIS{analysis}));
        %         selectivity_info.activity_info.(char(ANALYSIS{analysis})) = activity_info;
        selectivity_info.activity_info = activity_info;
        selectivity_info.num_iter.(char(ANALYSIS{analysis})) = NUM_ITER;
        
        fprintf([ANALYSIS{analysis}, ' ->selectivity_info saved.\n']);
        
        testName = 'Running Tests';
        multiWaitbar(testName, analysis/length(ANALYSIS), 'Color', [0.0 0.8 0.1] );
end


multiWaitbar(testName, 'Close');