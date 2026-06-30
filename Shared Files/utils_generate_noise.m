function [Mean_Data, Artificial_Var, Artificial_Rel_Var, Sample_Names] = utils_generate_noise(filename, n_replicates, noise_level)
% GENERATE_ARTIFICIAL_NOISE_DATA - Create artificial noisy data from the mean of experimental data and calculate its variance.

% Inputs:
%   filename     - Excel file with experimental data (e.g., 'ANDgate_data.xlsx')
%   n_replicates - The number of artificial noisy replicates to generate for each system.
%   noise_level  - The standard deviation of the Gaussian noise as a fraction of
%                  the mean value at each time point (e.g., 0.05 for 5% noise).

% Output:
%   Mean_Data        - Cell array containing the mean data trace for each unique sample.
%   Artificial_Var   - Cell array of the variance calculated from the artificial noisy data.
%   Artificial_Rel_Var - Cell array of the relative variance (variance / mean) from the artificial data.
%   Sample_Names     - Cell array of the unique sample names found in the file.

% Step 1: Read the data and calculate the mean of the real replicates
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
Data_table = readtable(filename, opts);

all_headers = Data_table.Properties.VariableNames;

% Find unique sample names by removing replicate suffixes like .1, .2, _1, etc.
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

% Initialize output cell arrays
Mean_Data = cell(1, length(Sample_Names));
Artificial_Var = cell(1, length(Sample_Names));
Artificial_Rel_Var = cell(1, length(Sample_Names));

min_threshold = 1e-6; % To prevent division by zero

% Step 2: Generate artificial data and calculate variance for each sample
for i = 1:length(Sample_Names)
    sample_name = Sample_Names{i};
    
    % Find all columns for the current sample to calculate the true mean
    replicate_mask = strcmp(base_names, sample_name);
    replicate_cols = valid_indices(replicate_mask);
    sample_data = table2array(Data_table(:, replicate_cols));
    
    % Calculate the mean trace from the real data
    mean_trace = mean(sample_data, 2, 'omitnan');
    Mean_Data{i} = mean_trace'; % Store as row vector
    
    % Create a matrix to hold the new noisy replicates
    artificial_replicates = zeros(length(mean_trace), n_replicates);
    
    for k = 1:n_replicates
        % Generate noise proportional to the mean signal at each time point
        % noise = sigma * N(0,1), where sigma = noise_level * signal_mean
        noise = (noise_level * mean_trace) .* randn(size(mean_trace));
        
        % Add noise to the mean trace to create one artificial replicate
        noisy_trace = mean_trace + noise;
        
        % Ensure non-negativity, as fluorescence can't be negative
        noisy_trace(noisy_trace < 0) = 0;
        
        artificial_replicates(:, k) = noisy_trace;
    end
    
    % Calculate variance across the new artificial replicates (dimension 2)
    if n_replicates > 1
        var_from_artificial = var(artificial_replicates, 0, 2);
    else
        % If only one replicate, variance is undefined. Use a small default.
        var_from_artificial = ones(size(mean_trace)) * min_threshold;
    end
    
    % Ensure variance is not zero to avoid issues in the objective function
    var_from_artificial(var_from_artificial < min_threshold) = min_threshold;
    
    Artificial_Var{i} = var_from_artificial'; % Store as row vector
    
    % Calculate relative variance (variance / mean)
    % Use a safe mean to avoid division by zero
    safe_mean = max(mean_trace, min_threshold);
    rel_var = var_from_artificial ./ (safe_mean.^2);
    
    Artificial_Rel_Var{i} = rel_var'; % Store as row vector
end

fprintf('Successfully generated %d artificial replicates with %.2f%% noise for %d systems.\n', ...
        n_replicates, noise_level*100, length(Sample_Names));

end
