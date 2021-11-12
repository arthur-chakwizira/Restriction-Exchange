%fit exponential to compartment population vs time and infer exchange rate
load('resex_mc_popn_1.mat', 'popn_1', 'popn_2', 'time')
% load('resex_mc_popn_1_gaussian.mat', 'popn_1', 'popn_2', 'time')

figure
plot(time, popn_1, time, popn_2)
N = popn_2(popn_1~=0);
t = time(popn_1~=0);
% N = popn_2(1:100);
% t = time(1:100);
logN = log(N);
p = polyfit(t', logN, 1);
disp(p)
tf = linspace(min(t), max(t), 1000);
Nf = polyval(p,tf);
figure
set(gcf, 'color', 'w');
set(gca,'box','on','linewidth',0.5,'layer','top', 'FontWeight', 'normal', 'FontSize', 12)
hold on
plot(tf, Nf, 'LineWidth', 2)
plot(t,logN,'--', 'LineWidth',2)
legend({'Fit', 'Data'})
xlabel('time [s]')
ylabel('Intracellular population (log)')
% end