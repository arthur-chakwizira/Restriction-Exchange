function signal = resex_mc_signal_general_km(m, xps)
%Returns signal according to full compartment model given protocol
b_values = xps.b;
q_vectors = xps.q; %this is now a cell array
n = xps.n;
dt =xps.dt;
%Model parameters
D_1 = m(1);
D_2 = m(2);
f_1 = m(3);
k = m(4);
%Derived parameters
f_2 = 1-f_1;
k_21= f_1*k;
k_12 = f_2*k;


K = [-k_12   k_21;   k_12   -k_21];
F = [f_1; f_2];


%predict signal using forward model

signal = zeros(n,1);

for c_n = 1:n
    qt = q_vectors{c_n};
    b = b_values(c_n);
    %% noticed that the scaling of q did not agree with the b-value
    qt = qt*sqrt(b/(sum(qt.^2)*dt));
    %%
    D = diag([D_1, D_2]);
    
    N_tp = numel(qt);
    exponential_product = 1;
    for i = 1:N_tp
        exponential_product = exponential_product*expm(K*dt - D.*(qt(i)^2)*dt);
    end
    
    signal(c_n) = ones(1,2)*exponential_product*F;
    
end
end
