function [Mean_Data, Std_Dev, Sample_Names, Time] = utils_get_mean_std(filename)
% GET_MEAN_AND_STD_FROM_RAW - Calculate mean and standard deviation from
% experimental data grouped by column headers.
% Input:
%   filename - Excel file with experimental data.
% Output:
%   Mean_Data    - Cell array of mean values for each unique sample.
%   Std_Dev      - Cell array of standard deviation for each sample.
%   Sample_Names - Cell array of unique sample names.
%   Time         - Time vector from the first column.

opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
Data_table = readtable(filename, opts);

Time = table2array(Data_table(:, 1));
all_headers = Data_table.Properties.VariableNames;

valid_headers = {};
valid_indices = [];
for i = 1:length(all_headers)
    header = all_headers{i};
    if ~contains(header, 'Unnamed', 'IgnoreCase', true) && ...
       ~contains(header, 'Time', 'IgnoreCase', true) && ...
       ~isempty(header)
        valid_headers{end+1} = header;
        valid_indices(end+1) = i;
    end
end

base_names = cellfun(@(h) regexprep(h, '[\._]\d+$', ''), valid_headers, 'UniformOutput', false);
Sample_Names = unique(base_names, 'stable');

Mean_Data = cell(1, length(Sample_Names));
Std_Dev = cell(1, length(Sample_Names));

for i = 1:length(Sample_Names)
    sample_name = Sample_Names{i};
    replicate_mask = strcmp(base_names, sample_name);
    replicate_cols = valid_indices(replicate_mask);
    sample_data = table2array(Data_table(:, replicate_cols));
    
    Mean_Data{i} = mean(sample_data, 2, 'omitnan')';
    
    if size(sample_data, 2) > 1
        Std_Dev{i} = std(sample_data, 0, 2, 'omitnan')';
    else
        Std_Dev{i} = zeros(1, size(sample_data, 1));
    end
end

end
