function [LeadsOut] = Neuro_DB_CheckLeads(curDir)
% Check_Leads

%----------------------------------------------------------------%
% Go to file location
%----------------------------------------------------------------%

cd(curDir);

LeadsOut = {};

%----------------------------------------------------------------%
% Load and Open Cheetah Log file
%----------------------------------------------------------------%
fid = fopen('CheetahLogFile.txt');

tline = fgetl(fid);
%----------------------------------------------------------------%
% Get Number of lines in the text file
%----------------------------------------------------------------%
count = 1;
while ischar(tline)
    tetout{count, 1} = tline;
    tline = fgetl(fid);
    count = count + 1;
end
%----------------------------------------------------------------%
% Close Cheetah Log file
%----------------------------------------------------------------%
fclose(fid);
%----------------------------------------------------------------%
% Create Directory location for NAS Neural data from mouse input
%----------------------------------------------------------------%
stTEt = strfind(tetout,'DISABLED');

disl = 0;
for line = 1:length(stTEt)
    if ~isempty(stTEt{line})
        disl = disl + 1;
    end
end

count = 1;
disStart = zeros(1,length(disl));
Line = zeros(1,length(disl));
for line = 1:length(stTEt)
    if ~isempty(stTEt{line})
        disStart(count) = stTEt{line};
        Line(count) = line;
        count = count + 1;
    end
end
%----------------------------------------------------------------%
% Create Cell array of disabled lines
%----------------------------------------------------------------%

for dc = 1:length(Line)
    tempLine = tetout{Line(dc)};
    tet = strfind(tempLine,'TT');
    LeadsOut{dc,1} = strcat('TT',tempLine(tet + 2));
    lead = strfind(tempLine,'has');
    LeadsOut{dc,2} = strcat('Lead_',tempLine(lead - 2));
end
    








