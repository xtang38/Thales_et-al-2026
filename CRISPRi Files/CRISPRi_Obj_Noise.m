function error = CRISPRi_Obj_Noise(p) 

global par Mean_Data Artificial_Var Artificial_Rel_Var Sample_Names 

% Number of conditions
num_conditions = length(Sample_Names);

tspan = 0:300:14400;  

options = odeset('RelTol',1e-5,'AbsTol',1e-7,'NonNegative',1:7);

%% CALCULATE ERROR FOR ALL CONDITIONS
error = 0;

for cond = 1:num_conditions
    
    % Initial conditions
    x0 = [0 0 0 35e-9 0 0 0 0 0 0];   
    
    % Set condition-specific trigger concentrations based on sample name
    condition_name = Sample_Names{cond};
     if strcmpi(condition_name, 'ON')
        par.Pcr = 0;    % crRNA promoter concentration (nM)
        par.Ptr = 0;    % tracrRNA promoter concentration (nM)
    
    elseif strcmpi(condition_name, 'C1')
        par.Pcr = 0.1*10^(-9);    % crRNA promoter concentration (nM)
        par.Ptr = 0.1*10^(-9);    % tracrRNA promoter concentration (nM)
    
    elseif strcmpi(condition_name, 'C2')
        par.Pcr = 0.25*10^(-9);    % crRNA promoter concentration (nM)
        par.Ptr = 0.25*10^(-9);    % tracrRNA promoter concentration (nM)
    
     end

     try
        % Solve ODE for this condition
        [t,x] = ode23s(@(t,x) CRISPRi_Model(t,x,p), tspan, x0, options);
        
        % Check for complex, NaN or Inf values in solution
        if ~isreal(x) || any(isnan(x(:))) || any(isinf(x(:)))
            error = 1e10;
            return;
        end

    % Extract mature GFP data (last state variable)
    GFP_sim = x(:,10)*10^6;  
    n_factor = 3.087043;
    GFP_sim = GFP_sim ./ n_factor;
    
    % Get experimental data for this condition

    Mean_exp = Mean_Data{cond};
    Var_w = Artificial_Rel_Var{cond};
    
    n_pts = min(length(GFP_sim), length(Mean_exp));
    g = GFP_sim(1:n_pts); g=g(:);
    m = Mean_exp(1:n_pts); m = m(:);
    v = Var_w(1:n_pts); v=v(:);
    error = error + sum(((g - m).^2)./v);

    catch ME
        % Bad parameter set — assign penalty and abort this evaluation
        fprintf('ODE failed for condition %s: %s\n', condition_name, ME.message);
        error = 1e10;
        return;
     end  
end