function [] = Neuro_DB_Mouse_Wrapper(M)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

MainLoc = 'G:\Tetrode_DATA\Days of Recording\';
mouseLoc = strcat(MainLoc,M,'\Neurophysiology');
NeuroDBLoc = strcat(MainLoc,'Neuron_Activity_Info_Database');

% Get Mouse Cells Previously analyzed
cd(NeuroDBLoc)
Ntempdir = cellstr(ls);
cellsProcessed = Ntempdir(3:end);
cellsProcDates = unique(cell2mat(cellfun(@(x) str2double(x(11:16)), cellsProcessed, 'UniformOutput', false)));


% Get list of mouse dates
cd(mouseLoc)
tempdir = cellstr(ls);
mouseDates = tempdir(3:end);
mouseAvailDates = cell2mat(cellfun(@(x) str2double(x), mouseDates, 'UniformOutput', false));


% Compare LIsts and get unAnalyzed list

dateNumindex = mouseAvailDates(~ismember(mouseAvailDates,cellsProcDates));
dateList = arrayfun(@(x) num2str(x), dateNumindex, 'UniformOutput', false);


% Compare LIsts and get unAnalyzed list

for di = 1:length(dateList)
    
    Neuro_DB_beta_v02(M,dateList{di});
    
end


end

