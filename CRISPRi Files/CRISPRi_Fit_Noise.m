clear all;
clc;
global par Mean_Data Artificial_Var Artificial_Rel_Var Sample_Names

filename = 'CRISPRdata.xlsx'; 
n_replicates = 6;               % Generate 6 artificial noisy replicates for each system
noise_level = 0.20;             % 40% Gaussian noise relative to the mean

par.Pytot = 0.5*10^(-9);  % Total target promoter concentration (nM)

[Mean_Data, Artificial_Var, Artificial_Rel_Var, Sample_Names] = utils_generate_noise(filename, n_replicates, noise_level);

% Parameter bounds
% Production rates (alpha)
Lalpha = 10^(-2);
Ualpha = 10;

% Binding rates (beta)
Lbind = 10^2;
Ubind = 10^6;

% Hybridization rates (gamma)
Lgamma = 10^3;
Ugamma = 10^7;

% Degradation rates (delta)
Ldeg = 10^(-4);
Udeg = 10^(-1);

% Promoter binding rate (omega)
Lomega = 10^3;
Uomega = 10^7;

% Ribosome binding/translation
Lrib = 10^(-4);
Urib = 10^(-2);

% Maturation rate
Lmat = 10^(-3);
Umat = 10^(-1);

% Lower bounds
lb = [Lalpha, Ldeg, Lgamma, Lalpha, Ldeg, Lgamma, Ldeg, ...
      Lomega, Lalpha, Ldeg, Lrib, Lrib, Lmat];

% Upper bounds  
ub = [Ualpha, Udeg, Ugamma, Ualpha, Udeg, Ugamma, Udeg, ...
      Uomega, Ualpha, Udeg, Urib, Urib, Umat];

p0 = (ub-lb).*rand(1,length(ub))+lb;

A = [];
b = [];
Aeq = [];
beq = [];
nlcon = [];

for iter = 1:30  % Number of iterations
    
    % GA options
    options = optimoptions('ga','PlotFcn','gaplotbestf','FunctionTolerance',1e-8, ...
                          'MaxGenerations', 100, 'PopulationSize',100);
    
    % Run genetic algorithm
    [p_opt,fval] = ga(@CRISPRi_Obj_Noise,length(p0),A,b,Aeq,beq,lb,ub,nlcon,options);
    
    % Calculate final error
    error = CRISPRi_Obj_Noise(p_opt);
    
    % Save results
    savefile = strcat('CRISPRonly_Noise20_RV_fitted',int2str(iter),'.mat');
    save(savefile,"p0","p_opt","error","fval")
    
    disp(['Iteration ' int2str(iter) ' completed. Error: ' num2str(error)]);
end
