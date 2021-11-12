function kappa_out = resex_mc_get_kappa_from_k(d, D, k_in, direction)
%Determine relationship between permeability and exchange rate
%Use non-linear least-squares curve fitting
%Using this model: k = ((d/a)*(1/kappa + d/(b*D))^-1
%From Tian X, et al. 2018
% close all

load('resex_mc_calib_v5.mat', 'calib')

if strcmp(direction, "from_1_to_2")
    kappa_k = calib.kappa_12.("d_" + num2str(d*1e6));
    if calib.D_1*1e9 ~= D*1e9; error('Input diffusivity does not match calibration'); end
elseif strcmp(direction, "from_2_to_1")
    kappa_k = calib.kappa_21.("d_" + num2str(d*1e6));
    if calib.D_2*1e9 ~= D*1e9; error('Input diffusivity does not match calibration'); end
else
    error('Invalid direction key. I do not understand what you want to do.')
end



%% Fit
kappa = kappa_k(:,1);
% k = kappa_k(:,2);
options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4,  'MaxIterations', 1e4,...
    'Display', 'off');
x0 = [6,10];
fun = @(x) ((d/x(1))*(1./kappa + d/(x(2)*D))).^-1 - kappa_k(:,2);
x = lsqnonlin(fun, x0, [], [], options);

% 
% figure
% set(gcf, 'color', 'w');
% set(gca,'box','on','linewidth',0.5,'layer','top', 'FontWeight', 'normal', 'FontSize', 12)
% 
% hold on
% plot(kappa,k, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7)
% 
% kappa_f = linspace(min(kappa), max(kappa));
% fun2 = @(x,kappa) ((d/x(1))*(1./kappa + d/(x(2)*D))).^-1;
% k_f = fun2(x, kappa_f);
% 
% plot(kappa_f, k_f, 'k-', 'LineWidth', 1)
% 
% legend({'Measurements', 'Fit'})
% % title('D = 1.2')
% % legend({'Measurements', 'Fit: D = 0.3e-9', 'Fit: D = 3e-9'})
% xlabel('Permeability [m/s]')
% ylabel('Exchange rate [/s]')

%get kappa given exchange
fun3 = @(k,x) (x(1)/(d*k) - d/(x(2)*D))^(-1);

kappa_out = fun3(k_in,x);
end
