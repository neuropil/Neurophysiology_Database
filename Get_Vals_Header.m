function [Values] = Get_Vals_Header(txtfile,searchStr)
%Get_Vals_Header: Gets info of interest from text file Header exported from
%neuralynx recording file
%   Detailed explanation goes here

ADBitVolts_str = cellfun(@(x) ~isempty(strfind(x,searchStr)), txtfile);

ADBitVolts_array = txtfile(ADBitVolts_str);

allstrs = strsplit(ADBitVolts_array{1}, ' ');
numstrs = allstrs(2:end);

numstrs(strcmp('',numstrs)) = [];

if length(numstrs) < 2
    outVal = numstrs{1};
    switch outVal
        case 'True'
            Values = 1;
        case 'False'
            Values = 0;
    end
else
    Values = cellfun(@(x) str2num(x), numstrs);
end

end

