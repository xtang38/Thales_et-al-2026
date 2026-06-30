clear all;
clc;
global par Mean_Data Var_Weight Sample_Names

filename = "ANDgate_data_Short.xlsx";

% Calculate mean and variance for each condition
[Mean_Data, Var_Weight, Sample_Names] = utils_get_mean_std(filename);

% p(1), p(2), p(3)  alpha      [nM·min⁻¹]
Lalpha= 1e-3;    Ualpha= 1;

% p(4)  beta_T   
Lhybrid= 1e-3;    Uhybrid= 10;

% p(6)  beta_leak 
Lbasal= 1e-4;    Ubasal= 1e-2;

%p(5), p(7), p(8) Hill Max rate
Lbleak= 1e-2;    Ubleak= 1;

% p(9) p(20) p(21) ki ke alpha_m   [min⁻¹]
Ltrans= 1e-4;    Utrans= 1;

% p(10)  p(11)   p(12) p(13) p(14) p(15) 
Ldeg= 1e-3;    Udeg= 1;

% p(16) n       
Lhill= 1;       Uhill= 4;

% p(17) p(18) p(19) K
LK= 1;       UK= 1e4; 

%p(22) scaling factor
Lscale= 1;      Uscale= 100;

lb = [Lalpha, Lalpha, Lalpha,  Lhybrid,  Lbleak,   Lbasal,  ...  % p1-p4
      Lbleak,   Lbleak,   Ltrans,  Ldeg,    ...  % p5-p8
      Ldeg,     Ldeg,     Ldeg,    Ldeg,    ...  % p9-p12
      Ldeg,     Lhill,    LK, LK, LK, Ltrans, Ltrans, Lscale];                   % p13-p15

ub = [Ualpha, Ualpha, Ualpha,  Uhybrid,  Ubleak,   Ubasal,  ...  % p1-p4
      Ubleak,   Ubleak,   Utrans,  Udeg,    ...  % p5-p8
      Udeg,     Udeg,     Udeg,    Udeg,    ...  % p9-p12
      Udeg,     Uhill,    UK,  UK, UK, Utrans, Utrans,  Uscale];                   % p13-p16

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
    [p_opt,fval] = ga(@ANDgate_Obj_Raw,length(p0),A,b,Aeq,beq,lb,ub,nlcon,options);
    
    % Calculate final error
    error = ANDgate_Obj_Raw(p_opt);
    
    % Save results
    savefile = strcat('ANDgate_Comp_fitted_',int2str(iter),'.mat');
    save(savefile,"p0","p_opt","error","fval")
    
    disp(['Iteration ' int2str(iter) ' completed. Error: ' num2str(error)]);
end
