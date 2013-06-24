function [cells2use_index] = Get_cellsIndex(cellNumber)

files = dir('*.mat');

cells2use = zeros(cellNumber,1);
cellcount = 1;
for stri = 1:length(files)
    if ~isempty(strfind(files(stri,1).name, 'clqual'))
        cells2use(cellcount,1) = stri;
        cellcount = cellcount + 1;
    end
end

list_quals = cell(length(cells2use),1);
cells2use_index = false(length(cells2use),1);
for qi = 1:length(list_quals)
    list_quals{qi} = files(cells2use(qi)).name;
    
    load(list_quals{qi});
    
    ID = CluSep.IsolationDist;
    LRat = CluSep.Lratio;
    
    if ID > 9 || LRat < 0.75 || isnan(LRat)
        cells2use_index(qi,1) = 1;
    else
        cells2use_index(qi,1) = 0;
    end
    
end



