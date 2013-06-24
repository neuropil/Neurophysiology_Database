function [Values] = Get_Vals_Header_regexp(txtfile,searchStr)
%Get_Vals_Header: Gets info of interest from text file Header exported from
%neuralynx recording file
%   Detailed explanation goes here

queries = {'-ADBitVolts','-InputRange','-ThreshVal','-InputInverted','-DualThresholding'};

lineIndex = cellfun(@(x) ~isempty(regexpi(x,searchStr)), txtfile, 'UniformOutput', true);

queryMatch = regexpi(txtfile{lineIndex},searchStr,'match');

queryIndex = find(strcmp(queryMatch{1},queries));

switch queryIndex
    case 1
        Values = str2double(regexp(txtfile{lineIndex},'0.\d*','match'));
    case 2 
        Values = str2double(regexp(txtfile{lineIndex},'\d*','match'));
    case 3
        Values = str2double(regexp(txtfile{lineIndex},'\d*','match'));
    case 4
        torf = regexp(txtfile{lineIndex},'[T|F]\w*','match');
        if strcmp(torf,'True')
            Values = 1;
        else
            Values = 0;
        end
    case 5
        torf = regexp(txtfile{lineIndex},'[T|F]\w*','match');
        if strcmp(torf,'True')
            Values = 1;
        else
            Values = 0;
        end
end









