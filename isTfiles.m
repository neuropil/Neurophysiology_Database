function [run_analysis] = isTfiles(curDir)
%ISTfiles Looks in directory for hte presence of .t files
%   Returns 0 if no .t files exist in directory and 1 if there are .t
%   files in the directory.  Absence of .t files indicates that the
%   recording day has not been clustered.

cd(curDir);

allfiles = dir('*.t');

if isempty(allfiles);
    run_analysis = 0;
else
    run_analysis = 1;
end

end

