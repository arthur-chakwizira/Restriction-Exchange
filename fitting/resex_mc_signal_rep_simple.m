function s = resex_mc_signal_rep_simple(m, xps)
% function s = resex_mc_signal_rep_simple(m, xps)
%
%
%

E_D  = m(1);
E_R  = m(2);
V = m(3);
k = m(4);

E = E_D + xps.Vomega*E_R;

s= exp(-xps.b.*E + (1/2)*V*(1-k*xps.Gamma).*xps.b.^2);
