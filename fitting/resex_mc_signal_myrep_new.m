function signal = resex_mc_signal_myrep_new(m, xps)
%Returns signal according to unified representation given protocol
%Does not use Gamma and Vomega, but the full expressions that these
%quantities approximate

%%
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
%Derived parameters
f_2 = 1-f_1;

R = d/2;

Qw = cell(size(q_vectors));
for c_n = 1:n
    q = q_vectors{c_n};
    [Qw{c_n}, ~] =  resex_mc_qt_to_qw(q, dt);
%     q = [q_vectors{c_n}  zeros(size(q_vectors{c_n}))];
%     Qw{c_n} = real(fftshift(fft(q)))*dt;
end



%Get diffusion spectrum
%Second-order approx
%         c = 7/1536;
%         Dw = w.^2*c*d^4/D_1;
%         full diffusion spectrum for cylinders


num_k = 200;
mu_k = my_besselzero(num_k);
Bk = 2*(R./mu_k).^2./(mu_k.^2 -1);
ak = (mu_k/R).^2;



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
        Dw(c_w) = sum(Bk.*ak.*D_1*ws^2./(ak.^2*D_1^2+ws^2), 'all');
    end
%     
%      Dw = w.^2*c*d^4/D_1;
    
    
    D1_t = (1/(pi*b))*trapz(w, Dw.*(qw));
    
    TE = (numel(qt)-1)*dt;
    t = (0:dt:TE);
    q4 = (1/b^2)*resex_mc__correlate(qt.^2, qt.^2, dt);
    h = 2*sum( exp(-k.*t).*q4 )*dt;
    
    
    E = f_1*D1_t + f_2*D_2;
    V = f_1*f_2*((D1_t - D_2)^2);
    signal(c_n)  = exp(-b*E + (1/2)*V*(h)*b^2);
    
end


end
