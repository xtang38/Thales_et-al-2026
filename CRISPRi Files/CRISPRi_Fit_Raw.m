
clear all;
clc;
global par Mean_Data Std_Dev Sample_Names

% Define inputs
filename = 'CRISPRdata_short.xlsx'; % Make sure this file is in the MATLAB path

par.Pytot = 0.5*10^(-9);  % Total target promoter concentration (nM)

% Calculate mean and variance for each condition from raw experimental data
[Mean_Data, Std_Dev, Sample_Names] = utils_get_mean_std(filename);

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

% Initialize random parameters
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
    [p_opt,fval] = ga(@CRISPRi_Obj_Raw,length(p0),A,b,Aeq,beq,lb,ub,nlcon,options);
    
    % Calculate final error
    error = CRISPRi_Obj_Raw(p_opt);
    
    % Save results
    savefile = strcat('CRISPRonly_Raw_noVar_Short_',int2str(iter),'.mat');
    save(savefile,"p0","p_opt","error","fval")
    
    disp(['Iteration ' int2str(iter) ' completed. Error: ' num2str(error)]);
end

