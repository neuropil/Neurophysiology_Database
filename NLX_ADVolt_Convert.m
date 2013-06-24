function [MicroVolts] = NLX_ADVolt_Convert(NTT_file)
%AD_Volt_Convert Converts AD Values for each lead on selected tetrode by
%ADBitValue stored in Header file of .ntt file.
%   INPUT: .ntt file (e.g.) 'TT1.ntt'
%        Example 1: Workspace variable : NTT_file = 'TT2.ntt'  
%                     Function call : AD_Volt_Convert(NTT_file)
%        Example 2: AD_Volt_Convert('TT2.ntt')  %% If .ntt file of interest 
%                   is located in Current Folder 
%   OUTPUT: Multidimension matrix : (m,n) equivalent to Samples(m,n)
%           with converted values
%
% DEPENDENCIES: Nlx2MatSpike.m

[~, ~, ~, ~, Samples, Header] = ...
            Nlx2MatSpike(NTT_file, [1 1 1 1 1], 1, 1, [] );

% Get ADBitVolts values from Header file       
ADBitVolts_str = Header(find(cellfun(@(x) ~isempty(strfind(x,'ADBitVolts')), Header)));
ADBitVolts_array = cell2mat(ADBitVolts_str);
[~, ADBit_Values] = strtok(ADBitVolts_array, ' ');
ADBitVolts = str2num(ADBit_Values);

convert_mat = repmat(ADBitVolts, [size(Samples, 1), 1, size(Samples, 3)]);
ConvertSamples = Samples .* convert_mat * 10^6;

MicroVolts = ConvertSamples;

end