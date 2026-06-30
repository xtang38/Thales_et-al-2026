function dx = CRISPRi_Model(t,x,p)

global par 

dx = zeros(10,1);

% State variables
crRNA = x(1);
tracrRNA = x(2);
gRNA = x(3);
dCAS9 = x(4);
Complex = x(5);
Py_bound = x(6);
M = x(7);
Mi = x(8);
G = x(9);
Gm = x(10);

% Free promoter 
Py_free = par.Pytot - Py_bound;

% Equation dcrRNA/dt
dx(1) = p(1)*par.Pcr - p(2)*crRNA - p(3)*crRNA*tracrRNA;

% Equation dtracrRNA/dt
dx(2) = p(4)*par.Ptr - p(5)*tracrRNA - p(3)*crRNA*tracrRNA;

% Equation dgRNA/dt
dx(3) =  p(3)*crRNA*tracrRNA - p(6)*gRNA*dCAS9 - p(7)*gRNA;

% Equation ddCAS9/dt
dx(4) = -p(6)*gRNA*dCAS9;

% Equation dComplex/dt
dx(5) = p(6)*gRNA*dCAS9 - p(8)*Complex*Py_free;

% Equation dPy_bound/dt
dx(6) = p(8)*Complex*Py_free;

% Equation dM/dt
dx(7) = p(9)*Py_free - p(10)*M - p(11)*M + p(12)*Mi;

% Equation dMi/dt
dx(8) = p(11)*M - p(12)*Mi;

% Equation dG/dt
dx(9) = p(12)*Mi - p(13)*G;

% Equation dGm/dt
dx(10) = p(13)*G;

return
