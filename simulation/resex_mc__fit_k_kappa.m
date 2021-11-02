function kappa_out = resex_mc__fit_k_kappa(d, D, k_in)
%Determine relationship between permeability and exchange rate
%Use non-linear least-squares curve fitting
%Using this model: k = ((d/a)*(1/kappa + d/(b*D))^-1
%From Tian X, et al. 2018
% close all

% d = 5; D = 1.2;



%% d =1, use dt = 1e-5;
% if d == 1 && D == 1.2
% D = 1.2e-9;
% d = d*1e-6;
% kappa_k=[
%     20e-6		20 ;
%     10e-6      3.6;
%     1.2e-6      1.3];
% end


%% d = 2, use dt = 1e-5
if d == 2 && D == 1.2
D = 1.2e-9;
R = 1e-6;
d = 2*R;
kappa_k=[
    20e-6		50 ;
    2.8e-6      7.2622;
    1.2e-6      2.6471];
end

%% d = 4, use dt = 1e-5
if d == 4 && D == 1.2
D = 1.2e-9;
R = 2e-6;
d = 2*R;
kappa_k=[
    20e-6		 26.1424;
    5.3553e-6    7.0608;
    2.4e-6     3];
end

%% d = 5, dt = 1e-4
if d == 5 && D == 1.2
   D = 1.2e-9;
   R = 2.5e-6;
   d = 2*R;
 %k_12
 kappa_k = [20e-6   20.3414;
     10e-6      9.6141;
     5e-6       4.8928
     2e-6       1.9842
     0.5e-6        0.5297;];
end


%% d = 6, dt = 1e-4
if d == 6 && D == 1.2
   D = 1.2e-9;
   R = 3e-6;
   d = 2*R;
 %k_12
 kappa_k = [20e-6   17;
     10e-6      8.7;
     5e-6       4.3;
     2e-6       1.7
     0.5e-6        0.43;];
end



%% d = 7, use dt = 1e-4
if d == 7 && D == 1.2
D = 1.2e-9;
R = 3.5e-6;
d = 2*R;
kappa_k=[40e-6      27.6;
    30e-6           23.1;
    20e-6	 12.4;
    15e-6       10.6;
    10e-6   7.19;
    5.0e-6   3.47 ;
    2e-6    1.4]; 
end

%% d = 8, use dt = 1e-4
if d == 8 && D == 1.2
D = 1.2e-9;
R = 4e-6;
d = 2*R;
kappa_k=[50e-6      31;
    30e-6           19.6;
    20e-6	 12.4;
    10e-6   6.19;
    5.0e-6   3. ;
    2e-6    1.27]; 
end


%% d = 9, dt = 1e-4
if d == 9 && D == 1.2
D = 1.2e-9;
R = 4.5e-6;
d = 2*R;
kappa_k=[50e-6         26.6;
    40e-6               22.8;
    30e-6                  17.4;
    20e-6			11.2 ;
    14e-6           7.81;
6e-6			3.4  ;
3e-6               1.53;
1.5e-6            0.76];
end


%% d = 10, dt = 1e-4
if d == 10 && D == 1.2
D = 1.2e-9;
R = 5e-6;
d = 2*R;
kappa_k=[50e-6         24.5;
    40e-6               20.2;
    30e-6                  14.1;
    20e-6			9.9288 ;
    14e-6           7.0226
6e-6			2.9292  ;
4e-6            2.0912 ;
3.5e-6			1.8039 ;
2e-6            0.9333 ;
0.35e-6			0.1772 ;
0.14e-6			0.0714  ;
0.07e-6			0.0315 ] ;
end

%% d = 12, use dt = 1e-4
if d ==12 && D == 1.2
D = 1.2e-9;
R = 6e-6;
d = 2*R;
kappa_k=[100e-6     33;
    50e-6           19.7;
    20e-6       8.05; %k_21
    10e-6      4.06;
    5e-6        2.06;
    2e-6        0.83;
    1e-6        0.4;
    ]; %k_12
end

%% d = 14, use dt = 1e-4
if d ==14 && D == 1.2
D = 1.2e-9;
R = 7e-6;
d = 2*R;
kappa_k=[100e-6     30.5;
    50e-6           17;
    20e-6       7; %k_21
    10e-6      3.6;
    5e-6        1.7;
    2e-6        0.75
    ]; %k_12
end


%% d = 16, dt = 1e-4
if d == 16 && D == 1.2
D = 1.2e-9;
R = 8e-6;
d = 2*R;
kappa_k=[
    100e-6      27.266;
    50e-6       15.0288;
    30e-6       8.915;
    20e-6       5.945;
    10e-6	 3.12; %k_12
    5e-6      1.73;
    2e-6        0.68]; %k_21
end

%% d = 18, dt = 1e-4
if d == 18 && D == 1.2
D = 1.2e-9;
R = 9e-6;
d = 2*R;
kappa_k=[100e-6   23.67;
    50e-6  13.4477;
    30e-6	 8.0377; %k_12
    20      5.4398;
    10      2.6388;
    5       1.4931;
    2       0.602]; %k_21
end


%% d = 20, dt = 1e-4
if d == 20 && D == 1.2
D = 1.2e-9;
R = 10e-6;
d = 2*R;
kappa_k=[100e-6     20;
    50e-6           12.1;
    20e-6           5.06;
    10e-6      2.48;
    5e-6            1.25;
    2e-6            0.5]; %k_12
end

%% Fit
kappa = kappa_k(:,1);
% k = kappa_k(:,2);
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
x0 = [6,10];
fun = @(x) ((d/x(1))*(1./kappa + d/(x(2)*D))).^-1 - kappa_k(:,2);
x = lsqnonlin(fun, x0, [], [], options);
% disp(x)

% kappa_f = linspace(min(kappa), max(kappa));
% fun2 = @(x,kappa) ((d/x(1))*(1./kappa + d/(x(2)*D))).^-1;
% k_f = fun2(x, kappa_f);

% figure
% set(gcf, 'color', 'w');
% set(gca,'box','on','linewidth',0.5,'layer','top', 'FontWeight', 'normal', 'FontSize', 12)
% 
% hold on
% plot(kappa,k, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7)
% plot(kappa_f, k_f, 'k-', 'LineWidth', 1)
% 
% legend({'Measurements', 'Fit'})
% % title('D = 1.2')
% % legend({'Measurements', 'Fit: D = 0.3e-9', 'Fit: D = 3e-9'})
% xlabel('Permeability [m/s]')
% ylabel('Exchange rate [/s]')
% 

%get kappa given exchange
fun3 = @(k,x) (x(1)/(d*k) - d/(x(2)*D))^(-1);

kappa_out = fun3(k_in,x);
end
