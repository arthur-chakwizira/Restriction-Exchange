%results collation simulations
% close all
% dock_figures(1)
%fwf
%Fixed d, varying k
%d = 10
true_k = [0   2     4      6       8        10];
k1 = [4        5    8.3     10.6    13.2   16.7];
k1_2 = [0 1.8 4.3 6.2 8.2 9.8];
d1 = [13.3  13.3    13.2    13.3    13.2     13]-3;



%Fixed k, varying d
%k = 10
true_d = [4       8     12       16     20];
d2 =      [3.3    10    15.7     20      115];
k2 =     [9.8     13    33       87     187 ];
d2_2 = [4.3 8 12.7 16 19.5];
k2_2 = [9.8 10.2 10.1 9.9 9.8];

if 1

subplot(3,2,3)
cla
set(gca, 'box', 'on')
set(gca, 'box', 'on', 'fontsize', 14)
hold on
plot(true_k, true_k, 'k-')
plot(true_k, k1_2, 's-', 'MarkerFaceColor', 'b', 'LineWidth', 2)
ylabel('Estimated k [s^{-1}]')
xlabel('Simulated k [s^{-1}]')


subplot(3,2,4)
cla
set(gca, 'box', 'on')
set(gca, 'box', 'on', 'fontsize', 14)
hold on
plot([0 true_d], [0 true_d], 'k-')
plot(true_d, d2_2, 's-', 'MarkerFaceColor', 'b', 'LineWidth', 2)
ylabel('Estimated d [\mum]')
xlabel('Simulated d [\mum]')

subplot(3,2,5)
cla
set(gca, 'box', 'on')
set(gca, 'box', 'on', 'fontsize', 14)
hold on
plot(true_k, d1, 's-', 'MarkerFaceColor', 'b', 'LineWidth', 2)
xlabel('Estimated k [s^{-1}]')
ylabel('Estimated d [\mum]')
ylim([0 11])



subplot(3,2,6)
cla
set(gca, 'box', 'on')
set(gca, 'box', 'on', 'fontsize', 14)
hold on
plot(true_d, k2_2, 's-', 'MarkerFaceColor', 'b', 'LineWidth', 2)
xlabel('Estimated d [\mum]')
ylabel('Estimated k [s^{-1}]')
ylim([0 11])


end




%vs b-value
%we use d4 k10
b = [2       5    10   15   20   25   30];
k_b = [9.9  10  9.8  9.9  10    10   10.1];
d_b = [4.2  3.3   3   2.9  2.9   3.2  4.4];
d_b2 =[4.1   3.9  4   3.9  3.9    4.1  4.];


%vs encoding time
%we use d4 k10
T = [  144    174    188    202    242  276  284  296];%ms
k_T = [10.4   10.6   10.7   10.9   9.6   12  11    11 ];
d_T = [4.4    4.2    1.7    4.0    4.6   4    4    0.3 ]; %?

k_T2 = [10.1 10. 10.2 9.9 9.6 10 10.1 9.9];
d_T2 = [4. 4.1 3.7 4. 4.2 4 4 3.9];

if 0

figure('Color', 'w')
subplot(1,3,2)
set(gca, 'box', 'on', 'fontsize', 20)
hold on
plot(b, ones(numel(b))*1, 'r-', 'LineWidth', 3, 'HandleVisibility', 'off')
plot(b, k_b/10, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
plot(b, d_b2/4, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8)
ylim([0 1.1])
xlabel('b [ms/\mu m^2]')
ylabel('Estimated/Simulated')
legend({'Exchange rate', 'Diameter'})

subplot(1,3,3)
set(gca, 'box', 'on', 'fontsize', 20)
hold on
plot(T, ones(numel(T), 1), 'r-', 'LineWidth', 3, 'HandleVisibility', 'off')
plot(T, k_T2/10, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
plot(T, d_T2/4, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8)
ylim([0 1.1])
xlabel('T [ms]')
ylabel('Estimated/Simulated')
legend({'Exchange rate', 'Diameter'})

end

%Figures to include:
%1. Discretisation of q-vector
%2. Substrate, protocol, gamma vomega
%3. Accuracy of k and d estimates (show fixed value on insert), contrast +
%fit goodness
%4. Variation of k and d with b, variation of k and d with T
%5. Precision
if 0
load('precision.mat')

pars(:,4) = pars(:,4)+0.5;

figure('Color', 'w')
subplot(1,2,1)
set(gca, 'box', 'on', 'fontsize', 16)
histogram(pars(:,4))
xlabel('k [s^{-1}]')
ylabel('Count')
hold on
plot([10 10], [0 25], 'r-', 'LineWidth', 3)
m_k = mean(pars(:,4));
std_k = std(pars(:,4));
legend({['\mu/\sigma ', ' = ', num2str(round(m_k, 1)), '/', num2str(round(std_k, 1))],...
    'Ground truth'})


subplot(1,2,2)
set(gca, 'box', 'on', 'fontsize', 16)
histogram(pars(:,5))
xlabel('d [\mum]')
ylabel('Count')
hold on
plot([4 4], [0 25], 'r-', 'LineWidth', 3)
m_d = mean(pars(:,5));
std_d = std(pars(:,5));

legend({['\mu/\sigma ', ' = ', num2str(round(m_d, 1)), '/', num2str(round(std_d, 1))],...
    'Ground truth'})
end
% subplot(3,2,3)
% set(gca, 'box', 'on', 'fontsize', 16)
% histogram(pars(:,1))
% xlabel('D_1 [\mum^2/ms]')
% ylabel('Count')
% 
% subplot(3,2,4)
% set(gca, 'box', 'on', 'fontsize', 16)
% histogram(pars(:,2))
% xlabel('D_2 [\mum^2/ms]')
% ylabel('Count')
% 
% subplot(3,2,5)
% set(gca, 'box', 'on', 'fontsize', 16)
% histogram(pars(:,3))
% xlabel('f_1')
% ylabel('Count')
