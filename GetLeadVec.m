function [LeadVec] = GetLeadVec(Tet2Use, TetLeadCells)
%GetLeadVec Summary of this function goes here
%   Detailed explanation goes here

tet = strcat('TT',num2str(Tet2Use));

tetrodes = TetLeadCells(:,1);
leads = TetLeadCells(:,2);

tetIndex = cellfun(@(x) strcmp(tet,x), tetrodes);

count = 1;
LeadVec = zeros(1,sum(tetIndex));
for ti = 1:length(tetIndex)
    if tetIndex(ti) == 1
        LeadVec(count) = str2double(leads{ti}(end)) + 1;
        count = count + 1;
    else
        continue
    end
end

