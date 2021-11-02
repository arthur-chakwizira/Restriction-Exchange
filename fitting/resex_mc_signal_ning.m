function signal = resex_mc_signal_ning(m, xps)

b = xps.b;
gammas = xps.Gamma;
%Model parameters
D_1 = m(1);
D_2 = m(2);
f_1 = m(3);
k = m(4);
%Derived parameters
f_2 = 1-f_1;

E = f_1*D_1 + f_2*D_2;
V = f_1*f_2*((D_1 - D_2).^2);
signal = exp(-b.*E + (1/2)*V.*(1-k*gammas).*b.^2);

end