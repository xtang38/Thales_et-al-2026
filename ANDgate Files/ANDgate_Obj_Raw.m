function error = ANDgate_Obj_Raw(p) 
% Fits all conditions simultaneously with one set of parameters

global par Mean_Data Var_Weight Sample_Names

% Number of conditions
num_conditions = length(Sample_Names);

%% SIMULATION SETUP
% Time span (adjust according to your experimental setup)
tspan = 0:10:240;  % 0 to 4 hours with 10 min intervals

% ODE solver options
options = odeset('RelTol',1e-5,'AbsTol',1e-7,'NonNegative',1:7);

% Target gene concentration (constant for all conditions)
par.Ctgt = 10;

%% CALCULATE ERROR FOR ALL CONDITIONS
error = 0;

for cond = 1:num_conditions
    
    % Initial conditions
    x0 = [0 0 0 0 0 0 0];  % [Ta Tb Tact Tgt Tgt_act G]
    
    % Set condition-specific trigger concentrations based on sample name
    condition_name = Sample_Names{cond};
     if strcmpi(condition_name, 'T1/T2')
        par.Cta = 20;
        par.Ctb = 20;
    
    elseif strcmpi(condition_name, 'T1')
        par.Cta = 20;
        par.Ctb = 0;
    
    elseif strcmpi(condition_name, 'T2')
        par.Cta = 0;
        par.Ctb = 20;
    
    elseif strcmpi(condition_name, 'Tar')
        par.Cta = 0;
        par.Ctb = 0;
    
     end

     try
        % Solve ODE for this condition
        [t,x] = ode23s(@(t,x) ANDgate_Model(t,x,p), tspan, x0, options);
        
        % Check for complex, NaN or Inf values in solution
        if ~isreal(x) || any(isnan(x(:))) || any(isinf(x(:)))
            error = 1e10;
            return;
        end
    
    % Extract mature GFP data (last state variable)
    GFP_sim = x(:,7);  % Mature GFP, convert to µM or appropriate units
    GFP_sim = GFP_sim.* p(22);
 
    n_factor = 5190.64;
    GFP_sim = GFP_sim ./ n_factor;
    
    % Get experimental data for this condition
    Mean_exp = Mean_Data{cond};
    Var_w = Var_Weight{cond};
    
    % Calculate error for this condition
    for i = 1:min(length(GFP_sim)-1, length(Mean_exp)-1)
        error = error + (GFP_sim(i+1) - Mean_exp(i))^2;
    end

    catch ME
        % Bad parameter set — assign penalty and abort this evaluation
        fprintf('ODE failed for condition %s: %s\n', condition_name, ME.message);
        error = 1e10;
        return;
     end  
end

