function [CellInfo] = Neuro_DB_beta(mouse,dateRec)

%%%%
% Incorporate Fieldtrip import and layout ideas
%%%%

BaseLoc = 'G:\Tetrode_DATA\Days of Recording\';
expectedNames = {'J303','J306','J311','J313','J314','J318','J319','J320'};

validateattributes(mouse,{'char'},{'nonempty'},mfilename,'mouseName', 1)
validateattributes(dateRec,{'char'},{'nonempty'},mfilename,'recdate', 2)

mouseName = validatestring(mouse,expectedNames,mfilename,'mouseName',1);
BehNLXLoc_all = strcat(BaseLoc,mouseName,'\Behavior_nlx\');
NeuroLoc_all = strcat(BaseLoc,mouseName,'\Neurophysiology\');

cd(BehNLXLoc_all);

poss_dates = cellstr(ls);

recdate = validatestring(dateRec,poss_dates,mfilename,'recdate', 2);

BehNLXLoc_date = strcat(BehNLXLoc_all,recdate);
NeuroLoc_date = strcat(NeuroLoc_all,recdate);

cd(NeuroLoc_date)
load('trialsTOincl.mat') % file name is trials_2_incl

run_analysis = isTfiles(NeuroLoc_date);

switch run_analysis
    
    case 0
        
        warndlg('Need to Cluster Data','User Error')
        
    case 1
        
        ftypes = {'t','ntt'};
        
        fileS = struct;
        for fi = 1:length(ftypes);
            
            fileLook = strcat('*.',ftypes{fi});
            dirFiles = dir(fileLook);
            
            for ni = 1:length(dirFiles)
                fileS.(ftypes{fi}){ni,1} = dirFiles(ni).name;
            end
            
        end
        
        % Run through leads
        
        for clustI = 1:length(fileS.t)
            
            % ADD MULTIBAR **** Examine function
            
            CellInfo = struct;
            
            getClustNum = fileS.t{clustI}(5);
            getTetNum = fileS.t{clustI}(3);
            
            cellName = strcat('Cell_',mouseName,'_',recdate,'_t',getTetNum,'_c',getClustNum);
            
            cd(NeuroLoc_date)
            
            NTT_file = char(strcat('TT',getTetNum,'.ntt'));
            
            load(strcat('TT',num2str(getTetNum),'cut.mat'))
            
            clustCut = eval(strcat('TT',num2str(getTetNum)));
            
            fprintf('CUT FILE for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            [~, ~, ~, ~, Samples, Header] = ...
                Nlx2MatSpike(NTT_file, [1 1 1 1 1], 1, 1, [] );
            
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
            
            disabledLeads = Neuro_DB_CheckLeads(NeuroLoc_date);
            
            % Remember Inactive leads will be 1:4 not 0:3 (NLX style)
            InactiveLeads = GetLeadVec(str2double(getTetNum),disabledLeads);
            
            cluSamples = Samples(:,:,clustCut == str2double(getClustNum));
            
            convert_mat = repmat(ADBitVolts, [size(cluSamples, 1), 1, size(cluSamples, 3)]);
            microVolts = cluSamples .* convert_mat * 10^6;
            
            if Inverted == 1;
                Clust_Waves = microVolts * -1;
            else
                Clust_Waves = microVolts;
            end
            
            load(strcat('TT',num2str(getTetNum),'_',getClustNum,'_clqual.mat'));
            
            ID = CluSep.IsolationDist;
            LR = CluSep.Lratio;
            
            allLeads = 1:4;
            disindex = allLeads(1:4 ~= InactiveLeads);
            
            fprintf('TIMESTAMPS for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            fprintf('DISABLED LEADS for cluster %d on tetrode %d in mouse %s recorded on %s evaluated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
                        
            features = struct;
            for fi = 1:length(disindex)
                
                tempWaves = squeeze(Clust_Waves(:,disindex(fi),:));
                
                features.Peak.(strcat('L',num2str(disindex(fi)))) =...
                    max(tempWaves)';
                
                features.Valley.(strcat('L',num2str(disindex(fi)))) =...
                    min(tempWaves)';
                
                features.Energy.(strcat('L',num2str(disindex(fi)))) =...
                    trapz(abs(tempWaves))';
                
                featsForPCA = horzcat(features.Peak.(strcat('L',num2str(disindex(fi)))),...
                    features.Valley.(strcat('L',num2str(disindex(fi)))),...
                    features.Energy.(strcat('L',num2str(disindex(fi)))));
                
                features.WavePC1.(strcat('L',num2str(disindex(fi)))) =...
                    pca(featsForPCA);
                
                %                 features.WaveSimIndex.(strcat('L',num2str(disindex(fi)))) =...
                %                     BrayCurtisIndex(tempWaves,disindex(fi));
                
                features.FSDE_Values.(strcat('L',num2str(disindex(fi)))) =...
                    FSDE_Method(tempWaves);
                
            end
            
            [features.WaveFitParams, features.WaveSumDS] = WaveFormFit(disindex,Clust_Waves);
            
            fprintf('WAVEFORM FEATURES for cluster %d on tetrode %d in mouse %s recorded on %s calculated...\n',...
                str2double(getClustNum), str2double(getTetNum), mouseName, recdate);
            
            % Experiment
            
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
            
            % Use date and mouse number to search Settings folder on
            % desktop; if file cannot be located search and copyfile from Z
            % drive: extract settings of interest
            localSettings = strcat('G:\All_Raw_Behavior\Settings\',mouse);
            ZSettings = strcat('Z:\Behavior\Behavior\Settings\John\',mouse);
            setName = strcat('settings_@Mix2afc_John_',mouse,'_',recdate,'a.mat');
            
            cd(localSettings)
            if ~exist(setName,'file')
                cd(ZSettings)
                if ~exist(mouseSet,'dir')
                    SLos = NaN;
                else
                    copyfile(setName,localSettings)
                    cd(localSettings)
                    load(setName)
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
                    SLos = mat2dataset(RngmrosS,'VarNames',{'Short','Long'});
                else
                    SLos = NaN;
                end
            end

            cd(BehNLXLoc_date)
            load(strcat('tb_',mouseName,'_',recdate,'.mat'))
            spkName = strcat('spk_tt',getTetNum,'_clust',getClustNum,'.mat');
            load(spkName)
            spktms = spk_fi/1000000;
            
            % For ISI calculation
            msSpktms = spktms*1000; % convert spike times to milliseconds
            spkIntervals = diff(msSpktms); 
            spkLogtimes = log(spkIntervals);
            % for plotting: hist(spkLogtimes,100);
            % for trouble shooting: ps = numel(find(spkIntervals < 1))/numel(spkIntervals)
            perISIviolate = numel(find(spkLogtimes < 0))/numel(spkLogtimes);
            
            % Taskbase stuff
            trStart = taskbase.start_nlx(1:end-1);
            
            epochNames = {'PreOdor','OdorSamp','Move','Reward','Baseline'};
            
            requiredTime = taskbase.req_time;
            
            midSLos = (SLos.Short + SLos.Long)/2;
            
            shortTrials = requiredTime < midSLos;
            longTrials = requiredTime > midSLos;
            
            preOEpoch = [taskbase.OdorPokeIn taskbase.DIO];
            odorSEpoch = [taskbase.DIO + 0.1 taskbase.OdorPokeOut];
            moveEpoch = [taskbase.OdorPokeOut taskbase.WaterPokeIn];
            rewEpoch = [taskbase.WaterPokeIn taskbase.WaterPokeIn + 0.75];
            allEpoch = [taskbase.OdorPokeIn taskbase.WaterPokeIn + 0.75];
            
            epochSS =  {preOEpoch,odorSEpoch,moveEpoch,rewEpoch,allEpoch};
            
            minStimID = min(taskbase.stimID);
            maxStimID = max(taskbase.stimID);
            stimIDs = minStimID:maxStimID;
            
            if stimIDs(1) == 8
                taskbase.stimID = taskbaase.stimID - 7;
            end
            
            behEvents = {'ovo','opo','wpi','wvo','wpo'};
            
            eventTimes = [taskbase.DIO taskbase.OdorPokeOut taskbase.WaterPokeIn taskbase.WaterDeliv taskbase.WaterPokeOut];
            %             sec_eventTimes = [taskbase.DIO taskbase.OdorPokeOut taskbase.WaterPokeIn taskbase.WaterDeliv taskbase.WaterPokeOut taskbase.NextOdorPokeIn];
            
            fr_DBs = {'EasyCorrect','HardCorrect','EasyError','HardError',...
                'LeftCorrect','RightCorrect','EasyAll','HardAll','LeftError',...
                'RightError','Fiftys'};
            
            TrialAnalyses = struct;
            TrialIDIndex = struct;
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
            
            % Cell Summary

            trialTypeNum = [0 1];            
            trials_to_in = trials_2_incl{clustI};
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
            
            % CHECK EPOCHS ANALYZED

            
            
            
            
            
            
            
            
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
            CellInfo.SpikeTimes = spktms;
            CellInfo.LRatio = LR;
            CellInfo.IsolationDistance = ID;
            CellInfo.Features = features;
            CellInfo.Experiment = Expermt;
            CellInfo.WaveIndex = WaveSimIndex;
            CellInfo.ShortLong_Info = SLos;
            CellInfo.ISIViolations = perISIviolate;
            
            cd('G:\Tetrode_DATA\Days of Recording\Neuron_Activity_Info_Database');
            
            save(cellName,'-struct','CellInfo');
            
        end
end



% Last things to do:
% Get short/long data in order RASTERS ##**!!**##
% Make sure all outputs in structure are accounted for
% Place fprintf points at each data analysis junction







% Future things
% 1. Calculate values for Autocorrelation
% 2. Turn CellSelectivity into switch construction




