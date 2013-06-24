function [Values] = stringSearch_beta(varargin)
%Get_Vals_Header: Gets info of interest from text file Header exported from
%neuralynx recording file
%
% REQUIRED Inputs : No labels, but order matters
% 'mouse': string variable: Example 'J303' 
% 'date': string variable: Example '111206' or 'Select' to bring up list
%
% OPTIONAL Inputs : REQUIRE labels preceding variables
% 'txtfile': string variable: 'CheetahLogFile' / 'Header' / 'Both'
%            DEFAULT: 'Both'
% 'query': cell array of strings or string variable: 'DISABLED',
%          'ADBitVolts','InputRange','InputInverted','ThreshVal',
%          'DualThresholding' or 'All' to search for all possible queries
%          DEFAULT: 'All'
% 'tetrode': integer variable: a vector or scalar of values 1:4
%            DEFAULT: 1:4
%
% EXAMPLES:
%
% Values = stringSearch_beta('J303','111210')
%     **     Output will derive all possible queries for all tetrodes
% Values = stringSearch_beta(...,'txtfile','Both')
%     **     Output will derive all possible queries for all tetrodes
% Values = stringSearch_beta(...,'txtfile','Both','query',{'DISABLED','InputRange'})
%     **     Output will derive selected queries for all possible tetrodes
% Values = stringSearch_beta(...,'txtfile','Both','query',{'DISABLED','InputRange'},'tetrode',2:3)
%     **     Output will derive selected queries for tetrodes 2 and 3


%==========================================================================
validQueries = {'DISABLED','ADBitVolts','InputRange','ThreshVal',...
    'InputInverted','DualThresholding'};
validTxtFiles = {'Header','CheetahLogFile','Both'};

p = inputParser;
p.addRequired('mouse',@(x) ischar(x));
p.addRequired('date', @(x) ischar(x));
p.addOptional('txtfile', 'Both', @(x) ismember(x,validTxtFiles));
p.addOptional('query', 'All', @(x) validAttributes(x));
p.addOptional('tetrode', 1:4, @(x) validAttributes(x))

% p.addParamValue('orientation', 'horiz', @(x) ismember(x, {'vert', 'horiz'}));
% p.addParamValue('transparency', 0.2, @(x) validateattributes(x, ...
% {'numeric'}, {'scalar', '>=', 0, '<=', 1}));

p.parse(varargin{:});

%-------------------------------------------------------------------------%
% Check if date is valid for mouse selected
%-------------------------------------------------------------------------%
BaseLoc = 'G:\Tetrode_DATA\Days of Recording\';
MouseLoc = strcat(BaseLoc,p.Results.mouse,'\Neurophysiology\');
cd(MouseLoc)

dir = cellstr(ls);
dateFiles = dir(3:end);

switch p.Results.date
    case 'Select'
        dateIndex = listdlg('PromptString','Select a Date',...
            'SelectionMode', 'single',...
            'ListSize', [100 250],...
            'ListString',dateFiles);
        
        date2use = dateFiles{dateIndex};
        
    otherwise
        
        if ~ismember(p.Results.date,dateFiles)
            disp('Invalid date selection: date for this mouse does not exist')
            dateIndex = listdlg('PromptSting','Select a Date',...
                'SelectionMode', 'single',...
                'ListString',dateFiles);
            
            date2use = dateFiles{dateIndex};
        else
            date2use = p.Results.date;
        end
end

%-------------------------------------------------------------------------%
% Go To Files Location
%-------------------------------------------------------------------------%
FileLoc = strcat(MouseLoc,date2use);
cd(FileLoc)

%-------------------------------------------------------------------------%
% Load Desired File
%-------------------------------------------------------------------------%

switch p.Results.txtfile
    case 'Header'
        HeaderAll = struct;
        for ti = 1:4
            temptet = strcat('TT',num2str(ti),'.ntt');
            [~, ~, ~, ~, ~, HeaderAll.(strcat('TT',num2str(ti)))] = ...
                Nlx2MatSpike(temptet, [1 1 1 1 1], 1, 1, [] );
        end
           
    case 'CheetahLogFile'
        fid = fopen('CheetahLogFile.txt');
        file = textscan(fid,'%s','Delimiter','\n');
        LeadData = file{1};
                
    case 'Both'
        HeaderAll = struct;
        for ti = 1:4
            temptet = strcat('TT',num2str(ti),'.ntt');
            [~, ~, ~, ~, ~, HeaderAll.(strcat('TT',num2str(ti)))] = ...
                Nlx2MatSpike(temptet, [1 1 1 1 1], 1, 1, [] );
        end
        fid = fopen('CheetahLogFile.txt');
        file = textscan(fid,'%s','Delimiter','\n');
        LeadData = file{1};
end

%-------------------------------------------------------------------------%
% Determine query vriable and convert if necessary
%-------------------------------------------------------------------------%
 
if ~iscell(p.Results.query)
    queries = p.Results.query;
    queries = cellstr(queries);
end

if any(strcmp(p.Results.query,'All')) && strcmp(p.Results.txtfile,'CheetahLogFile')
    queries = {'DISABLED'};
elseif any(strcmp(p.Results.query,'All')) && strcmp(p.Results.txtfile,'Both')
    queries = validQueries;
elseif strcmp(p.Results.txtfile,'Both') && any(strcmp(p.Results.query,'DISABLED'))
    queries = p.Results.query;
end

%-------------------------------------------------------------------------%
% Get analysis index based on query input
%-------------------------------------------------------------------------%

switch p.Results.txtfile
    case 'Header'
        analysis_flags = nan(5,1);
        for fi = 1:length(queries)
            analysis_flags(fi) = find(cellfun(@(x) strcmp(queries{fi},x),...
                validQueries));
        end
        
        analysis_flags = analysis_flags(~isnan(analysis_flags));
        
        validHeaderQs = 2:6;
        
        if ~all(ismember(analysis_flags, validHeaderQs))
            disp('Invalid query selection for Header');
            analysis_flags = input('Select values between range 2:6');
        end
        
    case 'CheetahLogFile'
        analysis_flags = find(cellfun(@(x) strcmp(queries,x),...
            validQueries));
        if analysis_flags ~= 1
            disp('Invalid Query in CheetahLogFile')
            analysis_flags = 1;
        end
        
    case 'Both' 
        analysis_flags = nan(6,1);
        for fi = 1:length(queries)
            analysis_flags(fi) = find(cellfun(@(x) strcmp(queries{fi},x),...
                validQueries));
        end
        
        analysis_flags = analysis_flags(~isnan(analysis_flags));
end

%-------------------------------------------------------------------------%
% Extract indices of interest for each query and tetrode
%-------------------------------------------------------------------------%

checkIndex = struct;
for qui = 1:length(analysis_flags)
    
    validCheck = 0;
    spaces_back = 1;
    
    if analysis_flags(qui) == 1;
        txtF = LeadData;
        switchCheck = 0;
    else
        switchCheck = 1;
    end
    
    switch switchCheck
        
        case 0
            
            while validCheck == 0
                
                firstTest = strcat(validQueries{analysis_flags(qui)}(1:spaces_back),'+\w*');
                
                lineIndex = cellfun(@(x) ~isempty(regexpi(x,firstTest)), txtF, 'UniformOutput', true);
                matchIndex = cellfun(@(x) regexpi(x,firstTest,'match'), txtF, 'UniformOutput', false);
                
                checkIndex.(validQueries{analysis_flags(qui)}) = find(lineIndex == 1);
                
                for li = 1:length(checkIndex.(validQueries{analysis_flags(qui)}))
                    checkline = matchIndex{checkIndex.(validQueries{analysis_flags(qui)})(li)}{1};
                    
                    if strcmp(checkline,validQueries{analysis_flags(qui)})
                        validCheck = 1;
                        continue
                    else
                        spaces_back = spaces_back + 1;
                        validCheck = 0;
                        break
                    end
                end
            end
            
        case 1
            
            for tti = 1:4
                txtF = HeaderAll.(strcat('TT',num2str(tti)));
                
                validCheck = 0;
                
                while validCheck == 0
                    
                    firstTest = strcat(validQueries{analysis_flags(qui)}(1:spaces_back),'+\w*');
                    
                    lineIndex = cellfun(@(x) ~isempty(regexpi(x,firstTest)), txtF, 'UniformOutput', true);
                    matchIndex = cellfun(@(x) regexpi(x,firstTest,'match'), txtF, 'UniformOutput', false);
                    
                    checkIndex.(validQueries{analysis_flags(qui)}).(strcat('TT',num2str(tti))) = find(lineIndex == 1);
                    
                    for li = 1:length(checkIndex.(validQueries{analysis_flags(qui)}).(strcat('TT',num2str(tti))))
                        checkline = matchIndex{checkIndex.(validQueries{analysis_flags(qui)}).(strcat('TT',num2str(tti)))(li)}{1};
                        
                        if strcmp(checkline,validQueries{analysis_flags(qui)})
                            validCheck = 1;
                            continue
                        else
                            spaces_back = spaces_back + 1;
                            validCheck = 0;
                            break
                        end
                    end
                end               
            end
    end  
end

%-------------------------------------------------------------------------%
% Return values of interest based on query input
%-------------------------------------------------------------------------%
tet2get = p.Results.tetrode;

for ayl2u = 1:length(analysis_flags)
    
    for ttFi = 1:length(tet2get)
        
        switch analysis_flags(ayl2u)
            
            case 1
                tempLines = checkIndex.(validQueries{analysis_flags(ayl2u)});
                tempTets = cell(length(tet2get),1);
                tempLeads = zeros(length(tet2get),1);
                for cl = 1:length(tempLines)
                    
                    temp_line = LeadData{tempLines(cl)};
                    
                    tt_ind = strfind(temp_line,'TT');
                    tt_out = temp_line(tt_ind:tt_ind + length('TT'));
                    tt = regexp(tt_out,'[0-9]','match');
                    
                    tempTets{cl,1} = strcat('TT',tt);
                    
                    chan_ind = strfind(temp_line,'Channel');
                    chan_out = temp_line(chan_ind(2):chan_ind(2) + length('Channel') + 1);
                    tempLeads(cl,1) = str2double(regexp(chan_out,'[0-9]','match'));
                end
 
                tetO = strcat('TT',num2str(tet2get(ttFi)));
                tIndex = cellfun(@(x) strcmp(tetO,x), tempTets);
                
                leadsOut = tempLeads(tIndex);
                
                Values.DisLeads.(strcat('TT',num2str(tet2get(ttFi)))) = leadsOut;

            case 2
                file2search = HeaderAll.(strcat('TT',num2str(tet2get(ttFi))));
                line2search = checkIndex.(validQueries{analysis_flags(ayl2u)}).(strcat('TT',num2str(tet2get(ttFi))));
                
                tempVals = str2double(regexp(file2search{line2search},'0.\d*','match'));
                Values.ADbitVolts.(strcat('TT',num2str(tet2get(ttFi)))) = tempVals(ttFi);
                
            case 3
                file2search = HeaderAll.(strcat('TT',num2str(tet2get(ttFi))));
                line2search = checkIndex.(validQueries{analysis_flags(ayl2u)}).(strcat('TT',num2str(tet2get(ttFi))));
                
                tempVals = str2double(regexp(file2search{line2search},'\d*','match'));
                Values.InputRange.(strcat('TT',num2str(tet2get(ttFi)))) = tempVals(ttFi);
                
            case 4
                file2search = HeaderAll.(strcat('TT',num2str(tet2get(ttFi))));
                line2search = checkIndex.(validQueries{analysis_flags(ayl2u)}).(strcat('TT',num2str(tet2get(ttFi))));
                
                tempVals = str2double(regexp(file2search{line2search}, '\d*','match'));
                Values.ThreshVals.(strcat('TT',num2str(tet2get(ttFi)))) = tempVals(ttFi);
                
            case 5
                file2search = HeaderAll.(strcat('TT',num2str(tet2get(ttFi))));
                line2search = checkIndex.(validQueries{analysis_flags(ayl2u)}).(strcat('TT',num2str(tet2get(ttFi))));
                
                torf = regexp(file2search{line2search},'[T|F]\w*','match');
                if strcmp(torf,'True')
                    Values.Inverted.(strcat('TT',num2str(tet2get(ttFi)))) = 1;
                else
                    Values.Inverted.(strcat('TT',num2str(tet2get(ttFi)))) = 0;
                end
                
            case 6
                file2search = HeaderAll.(strcat('TT',num2str(tet2get(ttFi))));
                line2search = checkIndex.(validQueries{analysis_flags(ayl2u)}).(strcat('TT',num2str(tet2get(ttFi))));
                
                torf = regexp(file2search{line2search},'[T|F]\w*','match');
                if strcmp(torf,'True')
                    Values.DualThreshold.(strcat('TT',num2str(tet2get(ttFi)))) = 1;
                else
                    Values.DualThreshold.(strcat('TT',num2str(tet2get(ttFi)))) = 0;
                end             
        end    
    end   
end

%-------------------------------------------------------------------------%
% Valid Attribute Function
%-------------------------------------------------------------------------%

function [valid] = validAttributes(inputvec)

validQueriesInput = {'DISABLED','ADBitVolts','InputRange','InputInverted','ThreshVal','DualThresholding','All'};

if ~iscell(inputvec) && ischar(inputvec)
    inputvec = cellstr(inputvec);
    test = 0;
elseif iscell(inputvec)
    test = 0;
end
     
if isnumeric(inputvec)
    test = 1;
end

switch test
    case 0
        test_valid = false(length(inputvec),1);
        for vi = 1:length(inputvec)
            test_valid(vi) = ismember(inputvec(vi),validQueriesInput);
        end
        
        if sum(test_valid) == length(test_valid)
            valid = 1;
        else
            valid = 0;
        end
    case 1
        if all(ismember(inputvec,1:4))
            valid = 1;
        else
            valid = 0;
        end
end

end



% Main Function End
end



