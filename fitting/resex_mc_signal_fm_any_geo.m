function signal = resex_mc_signal_fm_any_geo(m, xps)
%Returns signal according to full compartment model given protocol
%Valid for any geometry
b_values = xps.b;
q_vectors = xps.q; %this is now a cell array
n = xps.n;
dt =xps.dt;
%Model parameters
D_1 = m(1);
D_2 = m(2);
f_1 = m(3);
k = m(4);
d = m(5);
M = m(6);
N = m(7);
%Derived parameters
f_2 = 1-f_1;
k_21= f_1*k;
k_12 = f_2*k;


K = [-k_12   k_21;   k_12   -k_21];
F = [f_1; f_2];




Qw = cell(size(q_vectors)); %power spectrum
for c_n = 1:n
%     q = [q_vectors{c_n}  zeros(size(q_vectors{c_n}))];
%     Qw{c_n} = real(fftshift(fft(q)))*dt;
    
    qt = q_vectors{c_n}; 
    [Qw{c_n}, ~] = resex_mc_qt_to_qw(qt, dt); 
end



%Get diffusion spectrum
%Second-order approx
%         c = 7/1536;
%         Dw = w.^2*c*d^4/D_1;
%         full diffusion spectrum for cylinders


Bk = M*d^2;
ak = N/d^2;



%predict signal using forward model

signal = zeros(n,1);

for c_n = 1:n
    qw = Qw{c_n};
    qt = q_vectors{c_n};
    b = b_values(c_n);
    
    Nw = numel(qw);
    df = 1/(Nw*dt);
    if mod(Nw, 2) == 0; offset = 1; else; offset = 0; end
    w = 2*pi*df*(-floor(Nw/2):floor(Nw/2)-offset); %omega
    
    
    Dw = zeros(size(w));
    for c_w = 1:numel(w)
        ws = w(c_w);
        Dw(c_w) = Bk*ak*D_1*ws^2/(ak^2*D_1^2+ws^2);
    end
%     
%      Dw = w.^2*c*d^4/D_1;
    
    
%     D1_t = (1/(pi*b))*trapz(w, Dw.*(abs(qw).^2));
    D1_t = (1/(pi*b))*trapz(w, Dw.*(qw));
    D = diag([D1_t, D_2]);
    
    N_tp = numel(qt);
    exponential_product = 1;
    for i = 1:N_tp
        exponential_product = exponential_product*expm(K*dt - D.*(qt(i)^2)*dt);
    end
    
    signal(c_n) = ones(1,2)*exponential_product*F;
    
end
end
