function [CellInfo] = Neuro_DB_beta_v02(mouse,dateRec)

%%%%
% Incorporate Fieldtrip import and layout ideas
%%%%

% Location of all recording days for individual mice
BaseLoc = 'G:\Tetrode_DATA\Days of Recording\';
% Cell array of mice that have been recorded for experiments
expectedNames = {'J303','J306','J311','J313','J314','J318','J319','J320'};

% I don't remember how mfilename works in this construction
% Checks if mouseName variable is a character array
validateattributes(mouse,{'char'},{'nonempty'},mfilename,'mouseName', 1)
% Checks if recdate variable is a character array
validateattributes(dateRec,{'char'},{'nonempty'},mfilename,'recdate', 2)

% Checks whether selected mouseName is contained within Experimental mice
mouseName = validatestring(mouse,expectedNames,mfilename,'mouseName',1);

% If all input variables check out: then create folder indices
BehNLXLoc_all = strcat(BaseLoc,mouseName,'\Behavior_nlx\');
NeuroLoc_all = strcat(BaseLoc,mouseName,'\Neurophysiology\');

cd(BehNLXLoc_all);

% Get list of dates from mouse specified
poss_dates = cellstr(ls);

% Checks whether selected recDate is contained within mouse's date folder
recdate = validatestring(dateRec,poss_dates,mfilename,'recdate', 2);

% If all input variables check out: then navigate to specific Behavior_nlx
% folder and Neurophysiology folder (contains spike data and behavior)
BehNLXLoc_date = strcat(BehNLXLoc_all,recdate);
NeuroLoc_date = strcat(NeuroLoc_all,recdate);


% First navigate to Neurophysiology folder of Mouse on Recording day
cd(NeuroLoc_date)
% Load trials to include inforamtion (.mat file)
load('trialsTOincl.mat') % file name is trials_2_incl

% Check if tetrode files have been clustered
run_analysis = isTfiles(NeuroLoc_date);

% If the recording day has been clustered than run_analysis = 1
switch run_analysis
    
    case 0
        
        warndlg('Need to Cluster Data','User Error')
        
    case 1
        
        % File types to check for
        ftypes = {'t','ntt'};
        
        % Find all .t files and .ntt files
        % Each .t file indicates a unique cluster; 
        % Each .ntt file indicates a unique tetrode;
        fileS = struct;
        for fi = 1:length(ftypes);
            
            fileLook = strcat('*.',ftypes{fi});
            dirFiles = dir(fileLook); % creates a struct array of this file type
            
            % Sort file types into fileS struct
            for ni = 1:length(dirFiles)
                fileS.(ftypes{fi}){ni,1} = dirFiles(ni).name;
            end
            
        end
        
        totalNumclusters = length(fileS.t);
        numGoodqualClusters = length(trials_2_incl);
        [cells2use_index] = Get_cellsIndex(totalNumclusters);
        
        if numGoodqualClusters < sum(cells2use_index)
            difference = sum(cells2use_index) - numGoodqualClusters;
            for ii = 1:difference
                trials_2_incl{numGoodqualClusters + ii} = trials_2_incl{numGoodqualClusters};
            end
        end
        
        
        
        % Loop Through Clusters
        trials2useCount = 1;

        for clustI = 1:length(fileS.t)
            
            % ADD MULTIBAR **** Examine function
            
            % Initialize output structure
            
            

            
            
            
            
            
            
            
            
            

            
            CellInfo = struct;
            
            % Get Cluster number and Tetrode number from .t file name
            getClustNum = fileS.t{clustI}(5);
            getTetNum = fileS.t{clustI}(3);
            
            % Get Cell name from mouse,date,tetrode and cluster info
            cellName = strcat('Cell_',mouseName,'_',recdate,'_t',getTetNum,'_c',getClustNum);
            
            % Ensure that you're in Neurophysiology folder
            cd(NeuroLoc_date)
            
            % Create raw neurophysiology reference name
            NTT_file = char(strcat('TT',getTetNum,'.ntt'));
            
            % Load cluster index from cut file
            load(strcat('TT',num2str(getTetNum),'cut.mat'))
            
            % Create a usable variable in workspace from imported cut file
            clustCut = eval(strcat('TT',num2str(getTetNum)));
            
            fprintf('CUT FILE for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            % Extract all Samples data for waveforms and data from Header
            % for all Tetrodes and Leads
            [~, ~, ~, ~, Samples, Header] = ...
                Nlx2MatSpike(NTT_file, [1 1 1 1 1], 1, 1, [] );
            
            % Extract relevant data from Header
            ADBitVolts = Get_Vals_Header_regexp(Header,'-ADB[a-z]+');
            InputRange = Get_Vals_Header_regexp(Header,'-InputR[a-z]+');
            ThreshValues = Get_Vals_Header_regexp(Header,'-Thresh[a-z]+');
            Inverted = Get_Vals_Header_regexp(Header,'-InputI[a-z]+');
            DualThresh = Get_Vals_Header_regexp(Header,'-Dual[a-z]+');
            %             ADBitVolts = Get_Vals_Header(Header,'ADBitVolt');
            %             InputRange = Get_Vals_Header(Header,'InputRange');
            %             ThreshValues = Get_Vals_Header(Header,'ThreshVal');
            %             Inverted = Get_Vals_Header(Header,'InputInverted');
            %             DualThresh = Get_Vals_Header(Header,'DualThresholding');
            
            fprintf('HEADER for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            % Create variable to hold disabled lead information
            disabledLeads = Neuro_DB_CheckLeads(NeuroLoc_date);
            
            % Remember Inactive leads will be 1:4 not 0:3 (NLX style)
            InactiveLeads = GetLeadVec(str2double(getTetNum),disabledLeads);
            
            % Samples just for specific cluster
            cluSamples = Samples(:,:,clustCut == str2double(getClustNum));
            
            % Convert Samples into uVolts
            convert_mat = repmat(ADBitVolts, [size(cluSamples, 1), 1, size(cluSamples, 3)]);
            microVolts = cluSamples .* convert_mat * 10^6;
            
            % Check if Samples were inverted during recording and revert
            % the Samples if they were
            if Inverted == 1;
                Clust_Waves = microVolts * -1;
            else
                Clust_Waves = microVolts;
            end
            
            % Load cluster quality .mat file
            load(strcat('TT',num2str(getTetNum),'_',getClustNum,'_clqual.mat'));
            
            % Extract cluster quality data
            ID = CluSep.IsolationDist;
            LR = CluSep.Lratio;
            
            % Create index for disabled leads
            allLeads = 1:4;
            disindex = allLeads(1:4 ~= InactiveLeads);
            
            fprintf('TIMESTAMPS for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            fprintf('DISABLED LEADS for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            % Calculate features for recreation or rederivation of PCA
            features = struct;
            for fi = 1:length(disindex)
                
                % Reduce dimensions and extract Waveforms of interest
                tempWaves = squeeze(Clust_Waves(:,disindex(fi),:));
                
                % Waveform Structure
                WaveIndex.(strcat('L',num2str(disindex(fi)))) = tempWaves;
                
                % Peak
                features.Peak.(strcat('L',num2str(disindex(fi)))) =...
                    max(tempWaves)';
                % Valley
                features.Valley.(strcat('L',num2str(disindex(fi)))) =...
                    min(tempWaves)';
                % Energy
                features.Energy.(strcat('L',num2str(disindex(fi)))) =...
                    trapz(abs(tempWaves))';
                % Combine for WavePCA analysis
                featsForPCA = horzcat(features.Peak.(strcat('L',num2str(disindex(fi)))),...
                    features.Valley.(strcat('L',num2str(disindex(fi)))),...
                    features.Energy.(strcat('L',num2str(disindex(fi)))));
                % WavePC1
                features.WavePC1.(strcat('L',num2str(disindex(fi)))) =...
                    pca(featsForPCA);
                % Waveform Comparison (Bray-Curtis Index)
                % Keep off during dbug TAKES a while to compute 6/10/2013
                features.WaveSimIndex.(strcat('L',num2str(disindex(fi)))) =...
                    BrayCurtisIndex(tempWaves,disindex(fi));
                % First & Second Derivative analysis
                features.FSDE_Values.(strcat('L',num2str(disindex(fi)))) =...
                    FSDE_Method(tempWaves);
                
            end
            
            % Derive Gaussian fit parameters for positive and negative
            % component of waveform (Felsen and Thompson method)
            [features.WaveFitParams, features.WaveSumDS] = WaveFormFit(disindex,Clust_Waves);
            
            fprintf('WAVEFORM FEATURES for cluster %d on tetrode %d in mouse %s recorded on %s calculated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            % Determine Experiment index for mouse under examination           
            TF_2013 = {'J303','J306','J311','J313','J314'};
            TF_2014 = {'J318','J319','J320','J323','J330','J331'};
            TCF_2014 = {'J311','J313','J314','J318','J319'};
            
            expNames = {TF_2013,TF_2014,TCF_2014};
            expCategs = {'Thompson and Felsen 2013','Thompson and Felsen 2014',...
                'Thompson Costabile Felsen 2014'};
            
            expVec = false(1,3);
            for expi = 1:3
                expVec(expi) = ismember(mouse,(expNames{expi}));
            end
            
            Expermt = expCategs(expVec);

            % GET EVERYTHING DONE IN NEUROPHYSIOLOGY BEFORE MOVING ON
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%----------------------Behavior Section---------------------------%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Use date and mouse number to search Settings folder on
            % desktop; if file cannot be located search and copyfile from Z
            % drive: extract settings of interest
            localSettings = strcat('G:\All_Raw_Behavior\Settings\',mouse);
            ZSettings = strcat('Z:\Behavior\Behavior\Settings\John\',mouse);
            setName = strcat('settings_@Mix2afc_John_',mouse,'_',recdate,'a.mat');
            
            % Check if Short/Long trial were used
            cd(localSettings)
            % Check for settings file on computer drive
            if ~exist(setName,'file')
                % if it does not exist check the Z drive
                cd(ZSettings)
                if ~exist(setName,'file')
                    % if settings doesn't exist at ALL
                    % Set Short/Long trial Check to 0
                    SLos.Check = 0;
                else
                    % if the settings file does exist in the Z drive
                    copyfile(setName,localSettings)
                    cd(localSettings)
                    % Load settings file
                    load(setName)
                    setInfo = saved;
                    % Check and Extract Short Long time values
                    if strcmp(setInfo.OdorParameters_osp_mode,'Short/Long')
                        long_mros = setInfo.OdorParameters_m2;
                        long_rng = setInfo.OdorParameters_r2;
                        short_mros = setInfo.OdorParameters_m1;
                        short_rng = setInfo.OdorParameters_r1;
                        RngmrosS = [short_mros + short_rng , long_mros - long_rng];
                        SLos.Values = mat2dataset(RngmrosS,'VarNames',{'Short','Long'});
                        % If Short Long mode run then change check to 1
                        SLos.Check = 1;
                    else
                        % If not then change check to 0
                        SLos.Check = 0;
                    end
                end
            else
                load(setName)
                setInfo = saved;
                if strcmp(setInfo.OdorParameters_osp_mode,'Short/Long')
                    long_mros = setInfo.OdorParameters_m2;
                    long_rng = setInfo.OdorParameters_r2;
                    short_mros = setInfo.OdorParameters_m1;
                    short_rng = setInfo.OdorParameters_r1;
                    RngmrosS = [short_mros + short_rng , long_mros - long_rng];
                    SLos.Values = mat2dataset(RngmrosS,'VarNames',{'Short','Long'});
                    
                    SLos.Check = 1;
                else
                    SLos.Check = 0;
                end
            end
            
            % Go to behavior data location
            cd(BehNLXLoc_date)
            % Load behavior file
            load(strcat('tb_',mouseName,'_',recdate,'.mat'));
            spkName = strcat('spk_tt',getTetNum,'_clust',getClustNum,'.mat');
            % Load spike file
            load(spkName)
            % Transfer spike time data into seconds
            spktms = spk_fi/1000000;
            
            % For ISI calculation
            msSpktms = spktms*1000; % convert spike times to milliseconds
            spkIntervals = diff(msSpktms); 
            spkLogtimes = log(spkIntervals);
            % for plotting: hist(spkLogtimes,100);
            % for trouble shooting: ps = numel(find(spkIntervals < 1))/numel(spkIntervals)
            perISIviolate = numel(find(spkLogtimes < 0))/numel(spkLogtimes);
            
            ISIinfo.msSpktms = msSpktms;
            ISIinfo.spkIntervals = spkIntervals;
            ISIinfo.spkLogtimes = spkLogtimes;
            ISIinfo.ISIviolations = perISIviolate;
            
            % Taskbase stuff
            trStart = taskbase.start_nlx(1:end-1);
            
            % List of Epoch Names THESE Correspond to epochSS epochs
            epochNames = {'PreOdor','OdorSamp','Move','Reward','Baseline'};
            
            % Extract Required time variable for day of recording
            requiredTime = taskbase.req_time;
            
            % If Short Long mode used index short / long trials
            if SLos.Check
                midSLos = (SLos.Values.Short + SLos.Values.Long)/2;
                shortTrials = requiredTime < midSLos;
                longTrials = requiredTime > midSLos;
                % Short / Long trial types
                tType = {'All','Short','Long'};
                tIndicies = {~isnan(trStart),(shortTrials & ~isnan(trStart)),(longTrials & ~isnan(trStart))};
            end
            
            % Set up epochs used in behavior (timepoints to look between)
            preOEpoch = [taskbase.OdorPokeIn taskbase.DIO];
            odorSEpoch = [taskbase.DIO + 0.1 taskbase.OdorPokeOut];
            moveEpoch = [taskbase.OdorPokeOut taskbase.WaterPokeIn];
            rewEpoch = [taskbase.WaterPokeIn taskbase.WaterPokeIn + 0.75];
            allEpoch = [taskbase.OdorPokeIn taskbase.WaterPokeIn + 0.75];
            % Concatenate epochs
            epochSS =  {preOEpoch,odorSEpoch,moveEpoch,rewEpoch,allEpoch};
            
            % Determine if odor stimuli were taken from top or bottom box
            minStimID = min(taskbase.stimID);
            maxStimID = max(taskbase.stimID);
            stimIDs = minStimID:maxStimID;
            
            % Determine values used for stimID and ensure that they run
            % from 1 - 7
            if stimIDs(1) == 8
                taskbase.stimID = taskbaase.stimID - 7;
            end
            
            % Behavioral events
            % 'ovo': odor valve on
            % 'opo': odor poke out 
            % 'wpi': water poke in
            % 'wvo': water valve on 
            % 'wpo': water poke out
            behEvents = {'ovo','opo','wpi','wvo','wpo'};
            
            % Vectors for behavioral events extracted from taskbase; they
            % are times aligned to/with respect to trial start time
            eventTimes = [taskbase.DIO taskbase.OdorPokeOut taskbase.WaterPokeIn taskbase.WaterDeliv taskbase.WaterPokeOut];
            %             sec_eventTimes = [taskbase.DIO taskbase.OdorPokeOut taskbase.WaterPokeIn taskbase.WaterDeliv taskbase.WaterPokeOut taskbase.NextOdorPokeIn];
            
            % Odor / Outcome combination trial types
            % EasyCorrect: Rewarded trials for odorID = 95/5 80/20
            % HardCorrect: Rewarded trials for odorID = 60/40 
            % EasyError: Non-Rewarded trials for odorID = 95/5 80/20
            % HardError: Non-Rewarded trials for odorID = 60/40
            % LeftCorrect: Rewarded trials for odorID = 5 20 40
            % RightCorrect: Rewarded trials for odorID = 95 80 60
            % EasyAll: Rewarded and Non-R trials for odorID = 95/5 80/20
            % HardAll: Rewarded and Non-R trials for odorID = 60/40
            % LeftError: Non-Rewarded trials for odorID = 5 20 40
            % RightError: Non-Rewarded trials for odorID = 95 80 60
            % Fiftys: 50/50 trials

            fr_DBs = {'EasyCorrect','HardCorrect','EasyError','HardError',...
                'LeftCorrect','RightCorrect','EasyAll','HardAll','LeftError',...
                'RightError','Fiftys'};
            
            % Basic concept for Loop stucture of following code Loop for
            % each trial type and within Loop through trial sort type
            
            
            % Initialize Structures
            TrialAnalyses = struct; % Average Firing Rates
            TrialIDIndex = struct; % Trial Index
            % Check if there are short and long trials
            switch SLos.Check
                
                case 1 % There are short and long trials
                    % Loops through trial types (e.g. short/long) using structure names
                    for slm = 1:length(tType)
                        % Loops through trial sort types (e.g. LeftCorrect) using structure
                        % names
                        for dbi = 1:length(fr_DBs)
                            
                            switch dbi
                                case 1
                                    trialIndex = find(ismember(taskbase.stimID,[1 2 5 6]) &...
                                        taskbase.reward == 1 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 2
                                    trialIndex = find(ismember(taskbase.stimID,[3 4]) &...
                                        taskbase.reward == 1 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 3
                                    trialIndex = find(ismember(taskbase.stimID,[1 2 5 6]) &...
                                        taskbase.reward == 0 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 4
                                    trialIndex = find(ismember(taskbase.stimID,[3 4]) &...
                                        taskbase.reward == 0 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 5
                                    trialIndex = find(ismember(taskbase.stimID,[2 4 6]) &...
                                        taskbase.reward == 1 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 6
                                    trialIndex = find(ismember(taskbase.stimID,[1 3 5]) &...
                                        taskbase.reward == 1 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 7
                                    trialIndex = find(ismember(taskbase.stimID,[1 2 5 6]) &...
                                        tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 8
                                    trialIndex = find(ismember(taskbase.stimID,[3 4]) &...
                                        tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 9
                                    trialIndex = find(ismember(taskbase.stimID,[2 4 6]) &...
                                        taskbase.reward == 0 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 10
                                    trialIndex = find(ismember(taskbase.stimID,[1 3 5]) &...
                                        taskbase.reward == 0 & tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                case 11
                                    trialIndex = find(ismember(taskbase.stimID,7) &...
                                        tIndicies{slm});
                                    TrialIDIndex.(tType{slm}).(fr_DBs{dbi}) = trialIndex;
                                    
                            end
                            
                            % Calculate Raster values for each trial sort
                            % within each trial type ALIGNED TO ODOR POKE
                            % IN
                            [Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).opi, Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).TrialIndices,...
                                Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).RefTimes] = Neuro_DB_raster(spktms,...
                                trStart(trialIndex), taskbase.OdorPokeIn(trialIndex));
                            
                            % Get trial indices associated with spiketimes
                            % separated into blocks of different trials
                            % (NOT ACTUAL TRIAL THAT SPIKES OCCURRED IN)
                            tmp_trialNums = Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).TrialIndices;
                            
                            % REF_TRIALNUMS captures the ACTUAL TRIAL INDEX
                            % from where the SPIKETIME OCCURRED
                            ref_trialNums = zeros(length(tmp_trialNums),1);
                            for i = 1:length(tmp_trialNums)
                                ref_trialNums(i,1) = trialIndex(tmp_trialNums(i));
                            end
                            
                            % Save actual trial number where spikes
                            % occurred
                            Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).TrialIndices = ref_trialNums;
                            
                            % Ror each block of ACTUAL TRIAL INDICIES add
                            % REF_TIME (whichever event time used to derive time (starts with opi from above) and subtract ALIGN_TIME
                            % which is derived from: (trialstart time + event time(e.g opi))
                            
                            newReftimes = [];
                            for evi = 1:length(behEvents)
                                tempAllSpks = Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).opi;
                                tempReftimes = Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).RefTimes;
                                for trii = 1:length(trialIndex)
                                    refIndex = ref_trialNums == trialIndex(trii);
                                    tempSpikes = tempAllSpks(refIndex);
                                    zeroedSpikes = tempSpikes + tempReftimes(refIndex);
                                    alignSpikes = zeroedSpikes - (eventTimes(trialIndex(trii),evi) + trStart(trialIndex(trii)));
                                    newReftimes = [newReftimes; alignSpikes];
                                end
                                Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).(behEvents{evi}) = newReftimes;
                                [Ref_PSTH.(tType{slm}).(fr_DBs{dbi}).(behEvents{evi}),...
                                    PSTH_Info.(tType{slm}).(fr_DBs{dbi}).(behEvents{evi})] =...
                                    Neuro_DB_psth(newReftimes, Ref_SpkTimes.(tType{slm}).(fr_DBs{dbi}).TrialIndices);
                            end
                            
                            for epi = 1:length(epochSS)
                                
                                for tri = 1:length(trialIndex)
                                    tempTr = trialIndex(tri);
                                    stEp = epochSS{1,epi}(tempTr,1);
                                    enEp = epochSS{1,epi}(tempTr,2);
                                    epdur = enEp - stEp;
                                    
                                    if isempty(find(spktms > trStart(tempTr) + stEp & spktms < trStart(tempTr) + enEp))
                                        
                                        TrialAnalyses.(tType{slm}).(fr_DBs{dbi})(tri,epi) = NaN;
                                    else
                                        
                                        TrialAnalyses.(tType{slm}).(fr_DBs{dbi})(tri,epi) = numel(find(spktms >...
                                            trStart(tempTr) + stEp & spktms < trStart(tempTr) + enEp))/epdur;
                                    end
                                end
                            end
                            trial_DS = mat2dataset(TrialAnalyses.(tType{slm}).(fr_DBs{dbi}),'VarNames',epochNames);
                            TrialAnalyses.(tType{slm}).(fr_DBs{dbi}) = trial_DS;
                        end
                        
                    end
                    
                    
                case 0 % There are no short trials

                    for dbi = 1:length(fr_DBs)
                        
                        switch dbi
                            case 1
                                trialIndex = find(ismember(taskbase.stimID,[1 2 5 6]) &...
                                    taskbase.reward == 1 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 2
                                trialIndex = find(ismember(taskbase.stimID,[3 4]) &...
                                    taskbase.reward == 1 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 3
                                trialIndex = find(ismember(taskbase.stimID,[1 2 5 6]) &...
                                    taskbase.reward == 0 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 4
                                trialIndex = find(ismember(taskbase.stimID,[3 4]) &...
                                    taskbase.reward == 0 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 5
                                trialIndex = find(ismember(taskbase.stimID,[2 4 6]) &...
                                    taskbase.reward == 1 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 6
                                trialIndex = find(ismember(taskbase.stimID,[1 3 5]) &...
                                    taskbase.reward == 1 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 7
                                trialIndex = find(ismember(taskbase.stimID,[1 2 5 6]) &...
                                    ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 8
                                trialIndex = find(ismember(taskbase.stimID,[3 4]) &...
                                    ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 9
                                trialIndex = find(ismember(taskbase.stimID,[2 4 6]) &...
                                    taskbase.reward == 0 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 10
                                trialIndex = find(ismember(taskbase.stimID,[1 3 5]) &...
                                    taskbase.reward == 0 & ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                            case 11
                                trialIndex = find(ismember(taskbase.stimID,7) &...
                                    ~isnan(trStart));
                                TrialIDIndex.(fr_DBs{dbi}) = trialIndex;
                                
                        end
                        
                        [Ref_SpkTimes.(fr_DBs{dbi}).opi, Ref_SpkTimes.(fr_DBs{dbi}).TrialIndices,...
                            Ref_SpkTimes.(fr_DBs{dbi}).RefTimes] = Neuro_DB_raster(spktms,...
                            trStart(trialIndex), taskbase.OdorPokeIn(trialIndex));
                        
                        tmp_trialNums = Ref_SpkTimes.(fr_DBs{dbi}).TrialIndices;
                        
                        ref_trialNums = zeros(length(tmp_trialNums),1);
                        for i = 1:length(tmp_trialNums)
                            ref_trialNums(i,1) = trialIndex(tmp_trialNums(i));
                        end
                        
                        Ref_SpkTimes.(fr_DBs{dbi}).TrialIndices = ref_trialNums;
                        
                        % for each trial index add ref time and subtract align time
                        % which is derived from adding trialstart time and event
                        % time
                        
                        newReftimes = [];
                        for evi = 1:length(behEvents)
                            tempAllSpks = Ref_SpkTimes.(fr_DBs{dbi}).opi;
                            tempReftimes = Ref_SpkTimes.(fr_DBs{dbi}).RefTimes;
                            for trii = 1:length(trialIndex)
                                refIndex = ref_trialNums == trialIndex(trii);
                                tempSpikes = tempAllSpks(refIndex);
                                zeroedSpikes = tempSpikes + tempReftimes(refIndex);
                                alignSpikes = zeroedSpikes - (eventTimes(trialIndex(trii),evi) + trStart(trialIndex(trii)));
                                newReftimes = [newReftimes; alignSpikes];
                            end
                            Ref_SpkTimes.(fr_DBs{dbi}).(behEvents{evi}) = newReftimes;
                            [Ref_PSTH.(fr_DBs{dbi}).(behEvents{evi}),...
                                PSTH_Info.(fr_DBs{dbi}).(behEvents{evi})] =...
                                Neuro_DB_psth(newReftimes, Ref_SpkTimes.(fr_DBs{dbi}).TrialIndices);
                        end
                        
                        for epi = 1:length(epochSS)
                            
                            for tri = 1:length(trialIndex)
                                tempTr = trialIndex(tri);
                                stEp = epochSS{1,epi}(tempTr,1);
                                enEp = epochSS{1,epi}(tempTr,2);
                                epdur = enEp - stEp;
                                
                                if isempty(find(spktms > trStart(tempTr) + stEp & spktms < trStart(tempTr) + enEp))
                                    
                                    TrialAnalyses.(fr_DBs{dbi})(tri,epi) = NaN;
                                else
                                    
                                    TrialAnalyses.(fr_DBs{dbi})(tri,epi) = numel(find(spktms >...
                                        trStart(tempTr) + stEp & spktms < trStart(tempTr) + enEp))/epdur;
                                end
                            end
                        end
                        trial_DS = mat2dataset(TrialAnalyses.(fr_DBs{dbi}),'VarNames',epochNames);
                        TrialAnalyses.(fr_DBs{dbi}) = trial_DS;
                    end
                    
            end
            
            % Cell Summary

            trialTypeNum = [0 1];
            
            if cells2use_index(clustI)
                trials_to_in = trials_2_incl{trials2useCount};
                trials2useCount = trials2useCount + 1;
            else
                trials_to_in = 1:length(trStart);
            end
            
            
            align_ind = 3;
            window = [-3 3];
            spk_file = spk_fi;
            behav_file = taskbase;           
            
            for csi = 1:2
                
                switch csi
                    
                    case 1                                             
                        by_mixture_ratio_flag = trialTypeNum(csi);
                        
                        [Current.raster_info] = CellSummary_NDB(behav_file, spk_file, align_ind, window,...
                            by_mixture_ratio_flag, trials_to_in);
                        
                        Current.selectivity_info = Selectivity_Analysis_NDB(Current.raster_info);
                   
                    case 2
                        by_mixture_ratio_flag = trialTypeNum(csi);
                        
                        [Previous.raster_info] = CellSummary_NDB(behav_file, spk_file, align_ind, window,...
                            by_mixture_ratio_flag, trials_to_in);
                        
                        Previous.selectivity_info = Selectivity_Analysis_NDB(Previous.raster_info);
                end              
            end
            
            CellInfo.MouseName = mouseName;
            CellInfo.RecordDate = recdate;
            CellInfo.Tetrode = getTetNum;
            CellInfo.Cluster = getClustNum;
            CellInfo.ADBitVolts = ADBitVolts;
            CellInfo.InputRange = InputRange;
            CellInfo.ThreshValues = ThreshValues;
            CellInfo.Inverted = Inverted;
            CellInfo.DualThreshold = DualThresh;
            CellInfo.DisabledLeads = InactiveLeads;
            CellInfo.ClustWaves = Clust_Waves;
            CellInfo.ClusterIndex = clustCut;
            CellInfo.SpikeTimes = spktms;
            CellInfo.LRatio = LR;
            CellInfo.IsolationDistance = ID;
            CellInfo.Features = features;
            CellInfo.Experiment = Expermt;
            CellInfo.WaveIndex = WaveIndex;
            CellInfo.ShortLong_Info = SLos;
            CellInfo.ISIViolations = ISIinfo;
            CellInfo.CurrentTrials = Current;
            CellInfo.PreviousTrials = Previous;
            
            cd('G:\Tetrode_DATA\Days of Recording\Neuron_Activity_Info_Database');
            
            save(cellName,'-struct','CellInfo');
            
        end
end



% Last things to do:
% Make sure all outputs in structure are accounted for




% Future things
% 1. Calculate values for Autocorrelation
% 2. Turn CellSelectivity into switch construction