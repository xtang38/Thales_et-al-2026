function dx = ANDgate_Model(t,x,p)
global par
dx = zeros(7,1);

Ta      = max(x(1), 0);
Tb      = max(x(2), 0);
Tact    = max(x(3), 0);
Tgt_R   = max(x(4), 0);
Tgt_act = max(x(5), 0);
Mi      = max(x(6), 0);
G       = max(x(7), 0); 

% Parameters
% p(1)  = alpha:      transcription rate [nM·min⁻¹]
% p(2)  = alpha:      transcription rate [nM·min⁻¹]
% p(3)  = alpha:      transcription rate [nM·min⁻¹]

% p(4)  = beta_T:     RNA hybridization rate Ta+Tb->Tact [nM⁻¹·min⁻¹]
% p(5)  = beta_Tact:  hairpin opening efficiency of Tact (primary) [min⁻¹]
% p(6)  = beta_leak:  basal hairpin opening [min⁻¹]
% p(7)  = beta_Ta:    hairpin opening efficiency of Ta (leakage) [min⁻¹]
% p(8)  = beta_Tb:    hairpin opening efficiency of Tb (leakage) [min⁻¹]
% p(9)  = Ki:         Ribosome binding [min⁻¹]
% p(10)  = delta_G:   GFP degradation rate [min⁻¹]
% p(11)  = delta_Ta:  Ta RNA degradation rate [min⁻¹]
% p(12) = delta_Tb:   Tb RNA degradation rate [min⁻¹]
% p(13) = nu:         Tact dissociation rate [min⁻¹]
% p(14) = delta_TgtR: blocked mRNA degradation rate [min⁻¹]
% p(15) = delta_Tgta: open mRNA refolding/degradation rate [min⁻¹]
% p(16) = n:          Hill coefficient (shared) [-]
% p(17) = K:          Hill constant  [nM]
% p(18) = K:          Hill constant  [nM]
% p(19) = K:          Hill constant  [nM]
% p(20) = Ke:         Ribosome unbinding
% p(21) = K_tl:       GFP translation

H_Ta   = p(7)*(Ta^p(16))   / (Ta^p(16)   + p(17)^p(16));
H_Tb   = p(8)*(Tb^p(16))   / (Tb^p(16)   + p(18)^p(16));
H_Tact = p(5)*(Tact^p(16)) / (Tact^p(16) + p(19)^p(16));

opening_rate = Tgt_R*(H_Ta + H_Tb + H_Tact + p(6));

dx(1) = p(1)*par.Cta   + p(13)*Tact - p(4)*Ta*Tb - p(11)*Ta;
dx(2) = p(2)*par.Ctb   + p(13)*Tact - p(4)*Ta*Tb - p(12)*Tb;
dx(3) = p(4)*Ta*Tb     - p(13)*Tact;
dx(4) = p(3)*par.Ctgt - p(14)*Tgt_R - opening_rate;
dx(5) = opening_rate - p(15)*Tgt_act - p(9)*Tgt_act + p(20)*Mi;
dx(6) = p(9)*Tgt_act  - p(20)*Mi;
dx(7) = p(21)*Mi - p(10)*G;

return

