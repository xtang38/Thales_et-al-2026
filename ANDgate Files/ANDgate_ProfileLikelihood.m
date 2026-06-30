% profile_likelihood_ANDgate.m
% across all 9 optimisation sets (3 noise levels x 3 objective functions).
%   1. Loads the best parameter set (selected by re-evaluated SSE)
%   2. For each parameter p_i, fixes p_i at a series of values spanning
%      [lb_i, ub_i] and re-optimises all other parameters using fmincon
%   3. Plots the resulting PL curve and marks the 95% confidence threshold
%   4. Classifies each parameter as identifiable or non-identifiable
%   - One figure per optimisation set (9 figures total), each showing 21
%     PL curves in a grid layout
%   - Summary table comparing confidence interval widths across objective
%     functions
%   - Results saved as profile_likelihood_ANDgate_results.mat

clear; clc;
global par

%% ========================================================================
%  USER INPUTS — edit these to match your folder structure
%% ========================================================================

exp_data_file = 'ANDgate_Data_Short2.xlsx';

folder_names = {
    'Raw-noVariance',    'Raw-Variance',    'Raw-RV';
    'Noise-20-noVariance','Noise-20-Variance','Noise-20-RV';
    'Noise-40-noVariance','Noise-40-Variance','Noise-40-RV'
};
row_labels = {'Raw Data', 'Noise 20%', 'Noise 40%'};
col_labels  = {'No Variance', 'Variance', 'Relative Variance'};

% Parameter names (p22 scaling factor excluded)
param_names = {'\alpha_{Ta}', '\alpha_{Tb}', '\alpha_{Tgt}', ...
               '\beta_T', '\beta_{Tact}', '\beta_{leak}', ...
               '\beta_{Ta}', '\beta_{Tb}', 'K_i', '\delta_G', ...
               '\delta_{Ta}', '\delta_{Tb}', '\nu', '\delta_{TgtR}', ...
               '\delta_{Tgta}', 'n', 'K_1', 'K_2', 'K_3', ...
               'K_e', 'K_{tl}'};

% Parameter bounds (21 params, scaling factor p22 excluded)
lb = [1e-3, 1e-3, 1e-3,  1e-3,  1e-2,  1e-4, ...
      1e-2, 1e-2, 1e-4,  1e-3, ...
      1e-3, 1e-3, 1e-3,  1e-3, ...
      1e-3, 1,    1,     1,    1,    1e-4, 1e-4];

ub = [1,    1,    1,     10,   1,     1e-2, ...
      1,    1,    1,     1, ...
      1,    1,    1,     1, ...
      1,    4,    1e4,   1e4,  1e4,  1,    1];

n_profile_pts = 30;
delta_threshold = chi2inv(0.95, 1) / 2;

% ODE settings
tspan = 0:10:480;
ode_opts = odeset('RelTol',1e-4,'AbsTol',1e-6,'NonNegative',1:7,'MaxStep',10);
n_factor = 4497.26;

% Known constants
par.Ctgt = 10;

%% ========================================================================
%  LOAD EXPERIMENTAL DATA
%% ========================================================================
[exp_means, ~, sample_names, ~] = utils_get_mean_std(exp_data_file);
n_params = length(lb);

%% ========================================================================
%  MAIN LOOP
%% ========================================================================
results = struct();

for row = 1:3
    for col = 1:3
        folder = folder_names{row, col};
        set_label = sprintf('%s - %s', row_labels{row}, col_labels{col});
        fprintf('\n=== Processing: %s ===\n', set_label);

        mat_files = dir(fullfile(folder, '*.mat'));
        if isempty(mat_files)
            warning('No .mat files found in folder: %s', folder);
            continue;
        end

        best_sse = inf;
        best_p   = [];
        for fi = 1:length(mat_files)
            d = load(fullfile(folder, mat_files(fi).name));
            if isfield(d,'p_opt'), pv = d.p_opt(:)';
            elseif isfield(d,'p'), pv = d.p(:)';
            else, continue; end
            % Accept 21 or 22 parameter vectors
            if length(pv) == 22, pv = pv(1:21); end
            if length(pv) ~= n_params, continue; end
            sse_i = compute_sse_andgate(pv, sample_names, exp_means, tspan, ode_opts, n_factor);
            if sse_i < best_sse
                best_sse = sse_i;
                best_p   = pv;
            end
        end

        if isempty(best_p)
            warning('Could not load any valid parameter sets from: %s', folder);
            continue;
        end
        fprintf('  Best SSE = %.6f\n', best_sse);

        pl_values   = zeros(n_params, n_profile_pts);
        pl_grid     = zeros(n_params, n_profile_pts);
        ci_lower    = nan(1, n_params);
        ci_upper    = nan(1, n_params);
        identifiable = false(1, n_params);
        threshold = best_sse + delta_threshold;

        for pi = 1:n_params
            fprintf('  Parameter %d/%d: %s\n', pi, n_params, param_names{pi});

            if lb(pi) > 0
                grid_pts = logspace(log10(lb(pi)), log10(ub(pi)), n_profile_pts);
            else
                grid_pts = linspace(lb(pi), ub(pi), n_profile_pts);
            end
            pl_grid(pi,:) = grid_pts;

            fmin_opts = optimoptions('fmincon','Display','off', ...
                'Algorithm','interior-point', ...
                'MaxIterations',150, ...
                'MaxFunctionEvaluations',800, ...
                'OptimalityTolerance',1e-6, ...
                'StepTolerance',1e-8);

            for gi = 1:n_profile_pts
                fixed_val = grid_pts(gi);
                obj = @(p_free) compute_sse_andgate( ...
                    [p_free(1:pi-1), fixed_val, p_free(pi:end)], ...
                    sample_names, exp_means, tspan, ode_opts, n_factor);
                p0_free = [best_p(1:pi-1), best_p(pi+1:end)];
                lb_f    = [lb(1:pi-1), lb(pi+1:end)];
                ub_f    = [ub(1:pi-1), ub(pi+1:end)];
                try
                    [~, fval_pl] = fmincon(obj, p0_free, [], [], [], [], lb_f, ub_f, [], fmin_opts);
                    pl_values(pi, gi) = fval_pl;
                catch
                    pl_values(pi, gi) = 1e10;
                end
            end

            below = pl_values(pi,:) <= threshold;
            if any(below)
                idx_below = find(below);
                ci_lower(pi) = grid_pts(idx_below(1));
                ci_upper(pi) = grid_pts(idx_below(end));
                identifiable(pi) = (idx_below(1) > 1) && (idx_below(end) < n_profile_pts);
            end
        end

        results.(matlab.lang.makeValidName(folder)).best_p      = best_p;
        results.(matlab.lang.makeValidName(folder)).best_sse    = best_sse;
        results.(matlab.lang.makeValidName(folder)).pl_values   = pl_values;
        results.(matlab.lang.makeValidName(folder)).pl_grid     = pl_grid;
        results.(matlab.lang.makeValidName(folder)).ci_lower    = ci_lower;
        results.(matlab.lang.makeValidName(folder)).ci_upper    = ci_upper;
        results.(matlab.lang.makeValidName(folder)).identifiable = identifiable;
        results.(matlab.lang.makeValidName(folder)).threshold   = threshold;

        % Save partial results after each set in case of interruption
        save('profile_likelihood_ANDgate_results.mat', 'results', 'param_names', ...
             'folder_names', 'row_labels', 'col_labels', 'lb', 'ub');
        fprintf('  Partial results saved.\n');

        fig = figure('Position',[10,10,1600,1000],'Color','w');
        sgtitle(sprintf('Profile Likelihood — %s', set_label), ...
                'FontSize',13,'FontWeight','bold','Interpreter','tex');

        n_cols_pl = 6;
        n_rows_pl = ceil(n_params / n_cols_pl);

        for pi = 1:n_params
            ax = subplot(n_rows_pl, n_cols_pl, pi);
            hold(ax,'on');
            plot(ax, pl_grid(pi,:), pl_values(pi,:), 'b-', 'LineWidth', 1.8);
            yline(ax, threshold, 'r--', 'LineWidth', 1.2);
            [~, opt_idx] = min(pl_values(pi,:));
            plot(ax, pl_grid(pi, opt_idx), pl_values(pi, opt_idx), 'ko', ...
                 'MarkerFaceColor','k', 'MarkerSize', 5);
            set(ax, 'XScale','log', 'Box','off', 'FontSize', 8);
            xlabel(ax, param_names{pi}, 'Interpreter','tex', 'FontSize', 9);
            ylabel(ax, 'SSE', 'FontSize', 8);
            if identifiable(pi)
                title(ax, sprintf('CI: [%.2e, %.2e]', ci_lower(pi), ci_upper(pi)), ...
                      'FontSize', 7, 'Color', [0 0.5 0]);
            else
                title(ax, 'Non-identifiable', 'FontSize', 7, 'Color', [0.8 0 0]);
            end
        end

        saveas(fig, sprintf('PL_ANDgate_%s.png', matlab.lang.makeValidName(folder)));
        saveas(fig, sprintf('PL_ANDgate_%s.fig', matlab.lang.makeValidName(folder)));
    end
end

%% Summary table
fprintf('\n\n=== SUMMARY: Parameter Identifiability Across Objective Functions ===\n');
fprintf('%-20s | %-15s | %-15s | %-15s\n', 'Parameter', 'noVar', 'Var', 'RV');
fprintf('%s\n', repmat('-',1,70));
for pi = 1:n_params
    row_str = sprintf('%-20s', param_names{pi});
    for col = 1:3
        folder = folder_names{1, col};
        fn = matlab.lang.makeValidName(folder);
        if isfield(results, fn) && results.(fn).identifiable(pi)
            w = results.(fn).ci_upper(pi) - results.(fn).ci_lower(pi);
            row_str = [row_str, sprintf(' | CI=%.3e', w)];
        else
            row_str = [row_str, ' | Non-ident.   '];
        end
    end
    fprintf('%s\n', row_str);
end

save('profile_likelihood_ANDgate_results.mat', 'results', 'param_names', ...
     'folder_names', 'row_labels', 'col_labels', 'lb', 'ub');
fprintf('\nResults saved to profile_likelihood_ANDgate_results.mat\n');

%% ========================================================================
%  HELPER: compute SSE for AND gate
%% ========================================================================
function sse = compute_sse_andgate(p_vec, sample_names, exp_means, tspan, ode_opts, n_factor)
    global par
    sse = 0;
    % Use p(22) = 1 for SSE calculation (scaling factor excluded from PL)
    p_full = [p_vec(:)', 1];
    for ci = 1:length(sample_names)
        cname = sample_names{ci};
        if strcmpi(cname,'T1/T2')
            par.Cta = 20; par.Ctb = 20;
        elseif strcmpi(cname,'T1')
            par.Cta = 20; par.Ctb = 0;
        elseif strcmpi(cname,'T2')
            par.Cta = 0;  par.Ctb = 20;
        elseif strcmpi(cname,'Tar')
            par.Cta = 0;  par.Ctb = 0;
        end
        x0 = zeros(1,7);
        try
            [~,xsol] = ode23s(@(t,x) ANDgate_Model(t,x,p_full), tspan, x0, ode_opts);
            if ~isreal(xsol) || any(isnan(xsol(:))) || any(isinf(xsol(:)))
                sse = 1e10; return;
            end
            GFP_sim = xsol(:,7) * p_full(22) / n_factor;
            m = exp_means{ci};
            n_pts = min(length(GFP_sim), length(m));
            g = GFP_sim(1:n_pts); g = g(:);
            mv = m(1:n_pts); mv = mv(:);
            sse = sse + sum((g - mv).^2);
        catch
            sse = 1e10; return;
        end
    end
end
