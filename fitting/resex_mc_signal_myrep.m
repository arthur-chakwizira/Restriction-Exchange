function signal = resex_mc_signal_myrep(m, xps)
%Returns signal according to unified representation given protocol
b = xps.b;
gammas = xps.Gamma;
vomegas = xps.Vomega;
%Model parameters
D_1 = m(1);
D_2 = m(2);
f_1 = m(3);
k = m(4);
d = m(5);
%Derived parameters
f_2 = 1-f_1;
c = 7/1536;
R = c*d^4/D_1; %restriction coefficient

E = f_1*R.*vomegas + f_2*D_2;
V = f_1*f_2*((R.*vomegas - D_2).^2);
signal = exp(-b.*E + (1/2)*V.*(1-k*gammas).*b.^2);
end
