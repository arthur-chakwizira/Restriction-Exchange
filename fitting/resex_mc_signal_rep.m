function signal = resex_mc_signal_rep(m, xps)
%Returns signal according to unified representation given protocol
b = xps.b;
gammas = xps.Gamma;
vomegas = xps.Vomega;
%Model parameters

E_D  = m(1);
E_R  = m(2);
V_D = m(3);
C_DR = m(4);
V_R = m(5);
k = m(6);

E = E_D + vomegas*E_R;
V = V_D + 2*vomegas*C_DR + vomegas.^2*V_R;
signal = exp(-b.*E + (1/2)*V.*(1-k*gammas).*b.^2);
end
