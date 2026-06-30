% profile_likelihood_CRISPRonly.m
% model across all 9 optimisation sets (3 noise levels x 3 objective functions).
%   1. Loads the best parameter set (selected by re-evaluated SSE, same
%      logic as in plot_CRISPRonly_GA_results.m)
%   2. For each parameter p_i, fixes p_i at a series of values spanning
%      [lb_i, ub_i] and re-optimises all other parameters using fmincon
%   3. Plots the resulting PL curve and marks the 95% confidence threshold
%   4. Classifies the parameter as identifiable or non-identifiable based
%      on whether the profile crosses the threshold on both sides
%   - One figure per optimisation set (9 figures total), each showing 13
%     PL curves in a grid layout
%   - A summary table comparing confidence interval widths across the 3
%     objective functions for each noise level
%   - Results saved as profile_likelihood_CRISPRonly_results.mat
%   - CRISPRi_Model.m and utils_get_mean_std.m in the same folder
%   - The 9 results folders with 30 .mat files each in the same folder
%   - Optimization Toolbox (for fmincon)

clear; clc;
global par

%% ========================================================================
%  USER INPUTS — edit these to match your folder structure
%% ========================================================================

% Path to experimental data file
exp_data_file = 'CRISPRdata.xlsx';

% Results folders: rows = noise levels, cols = objective functions
folder_names = {
    'Raw-noVar',     'Raw-Var',     'Raw-RV';
    'Noise-20-noVar','Noise-20-Var','Noise-20-RV';
    'Noise-40-noVar','Noise-40-Var','Noise-40-RV'
};
row_labels = {'Raw Data', 'Noise 20%', 'Noise 40%'};
col_labels  = {'No Variance', 'Variance', 'Relative Variance'};

% Parameter names (for plot labels)
param_names = {'\alpha_{cr}', '\delta_{cr}', '\gamma_{hyb}', ...
               '\alpha_{tr}', '\delta_{tr}', '\gamma_{gRNA}', ...
               '\delta_{gRNA}', '\omega', '\alpha_m', '\delta_m', ...
               'k_{rib}', 'k_e', 'k_{mat}'};

% Parameter bounds — matched to GA fitting script
% p(1)=alpha_cr, p(2)=delta_cr, p(3)=gamma_hyb, p(4)=alpha_tr, p(5)=delta_tr
% p(6)=gamma_gRNA, p(7)=delta_gRNA, p(8)=omega, p(9)=alpha_m, p(10)=delta_m
% p(11)=k_rib, p(12)=k_e, p(13)=k_mat
lb = [1e-2,  1e-4,  1e3,  1e-2,  1e-4,  1e3,  1e-4,  1e3,  1e-2,  1e-4,  1e-4,  1e-4,  1e-3];
ub = [10,    1e-1,  1e7,  10,    1e-1,  1e7,  1e-1,  1e7,  10,    1e-1,  1e-2,  1e-2,  1e-1];

% Number of profile points per parameter (more = smoother curve, slower)
n_profile_pts = 30;

% Confidence threshold: chi2 critical value for 1 DOF at 95%
% Delta_chi2 = chi2inv(0.95,1)/2 = 1.9208 added to minimum SSE
delta_threshold = chi2inv(0.95, 1) / 2;

% ODE settings
tspan = 0:300:14400;
ode_opts = odeset('RelTol',1e-5,'AbsTol',1e-7,'NonNegative',1:10);
n_factor = 3.087043;

% Known model constants
par.Pytot = 0.5e-9;

%% ========================================================================
%  LOAD EXPERIMENTAL DATA
%% ========================================================================
[exp_means, ~, sample_names, ~] = utils_get_mean_std(exp_data_file);
n_params = length(lb);

%% ========================================================================
%  MAIN LOOP: iterate over all 9 sets
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
            if length(pv) ~= n_params, continue; end
            sse_i = compute_sse_crispr(pv, sample_names, exp_means, tspan, ode_opts, n_factor);
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

        pl_values    = zeros(n_params, n_profile_pts);
        pl_grid      = zeros(n_params, n_profile_pts);
        ci_lower     = nan(1, n_params);
        ci_upper     = nan(1, n_params);
        identifiable = false(1, n_params);

        threshold = best_sse + delta_threshold;

        fmin_opts = optimoptions('fmincon', 'Display','off', ...
            'Algorithm','interior-point', ...
            'MaxFunctionEvaluations', 5000, ...
            'OptimalityTolerance', 1e-8);

        for pi = 1:n_params
            fprintf('  Parameter %d/%d: %s\n', pi, n_params, param_names{pi});

            % Log-spaced grid from lb to ub
            grid_pts = logspace(log10(lb(pi)), log10(ub(pi)), n_profile_pts);
            pl_grid(pi,:) = grid_pts;

            for gi = 1:n_profile_pts
                fixed_val = grid_pts(gi);

                % Objective: SSE with p_i fixed at fixed_val
                obj = @(p_free) compute_sse_crispr( ...
                    [p_free(1:pi-1), fixed_val, p_free(pi:end)], ...
                    sample_names, exp_means, tspan, ode_opts, n_factor);

                % Initial guess and bounds with p_i removed
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

            % Determine confidence interval
            below = pl_values(pi,:) <= threshold;
            if any(below)
                idx_below = find(below);
                ci_lower(pi) = grid_pts(idx_below(1));
                ci_upper(pi) = grid_pts(idx_below(end));
                % Identifiable if profile rises above threshold on both sides
                identifiable(pi) = (idx_below(1) > 1) && (idx_below(end) < n_profile_pts);
            end
        end

        fn = matlab.lang.makeValidName(folder);
        results.(fn).best_p       = best_p;
        results.(fn).best_sse     = best_sse;
        results.(fn).pl_values    = pl_values;
        results.(fn).pl_grid      = pl_grid;
        results.(fn).ci_lower     = ci_lower;
        results.(fn).ci_upper     = ci_upper;
        results.(fn).identifiable = identifiable;
        results.(fn).threshold    = threshold;

        fig = figure('Position',[10,10,1400,900],'Color','w');
        sgtitle(sprintf('Profile Likelihood — %s', set_label), ...
                'FontSize',13,'FontWeight','bold','Interpreter','tex');

        n_cols_pl = 5;
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

        saveas(fig, sprintf('PL_CRISPRonly_%s.png', matlab.lang.makeValidName(folder)));
        saveas(fig, sprintf('PL_CRISPRonly_%s.fig', matlab.lang.makeValidName(folder)));
        close(fig);
    end
end

%% ========================================================================
%  SUMMARY TABLE: CI widths across objective functions
%% ========================================================================
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

save('profile_likelihood_CRISPRonly_results.mat', 'results', 'param_names', ...
     'folder_names', 'row_labels', 'col_labels', 'lb', 'ub');
fprintf('\nResults saved to profile_likelihood_CRISPRonly_results.mat\n');

%% ========================================================================
%  LOCAL FUNCTION — must be at the end of the script
%% ========================================================================
function sse = compute_sse_crispr(p_vec, sample_names, exp_means, tspan, ode_opts, n_factor)
    global par
    sse = 0;
    for ci = 1:length(sample_names)
        cname = sample_names{ci};
        if strcmpi(cname,'ON')
            par.Pcr = 0; par.Ptr = 0;
        elseif strcmpi(cname,'C1')
            par.Pcr = 0.1e-9; par.Ptr = 0.1e-9;
        elseif strcmpi(cname,'C2')
            par.Pcr = 0.25e-9; par.Ptr = 0.25e-9;
        end
        x0 = [0 0 0 35e-9 0 0 0 0 0 0];
        try
            [~,xsol] = ode23s(@(t,x) CRISPRi_Model(t,x,p_vec), tspan, x0, ode_opts);
            if ~isreal(xsol) || any(isnan(xsol(:))) || any(isinf(xsol(:)))
                sse = 1e10; return;
            end
            GFP_sim = xsol(:,10)*1e6 / n_factor;
            m  = exp_means{ci};
            n_pts = min(length(GFP_sim), length(m));
            g  = GFP_sim(1:n_pts); g = g(:);
            mv = m(1:n_pts);       mv = mv(:);
            sse = sse + sum((g - mv).^2);
        catch
            sse = 1e10; return;
        end
    end
end
