function resex_mc_study_k_kappa()
close all
d = 5e-6;%linspace(5e-6, 15e-6); %diameter
a = 6; %constant
b = 6e3; %constant
kappa =5e-6; %permeability in m/s


D = linspace(0, 2)*1e-9;
k = ((d/a).*(1/kappa + d./(b.*D))).^-1;

figure('Color', 'w')
set(gca,'box','on','linewidth',1, 'FontWeight', 'normal', 'FontSize', 20)
grid minor
hold on
plot(D*1e9, k, 'k-', 'LineWidth', 3)

ylabel('$k$ [s$^{-1}$]', 'interpreter', 'latex', 'fontsize', 30)
    xlabel('$D$ [$\mu ^2$/ms]', 'interpreter', 'latex', 'fontsize', 30)

end