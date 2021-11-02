function resex_mc_plot_s_vs_gamma_vs_vomega(S, xps, f)
% f = figure('Color', 'w');
% h = axes(f);
% hold(h, 'on')
%extract signals generated with same waveform and display them
[gam_vom, groups, ind] = unique([round(xps.Gamma*1e3, 4) round(1e3./sqrt(xps.Vomega), 4)], 'rows');
sigs = NaN(numel(groups), 20); %assuming at most 10 b-values

for c_g = 1:numel(groups)
    tmp_s = S(ind == (c_g)); %signals with same waveform
    tmp_b = xps.b(ind == (c_g));
    [b_s, b_groups, b_ind] = unique(tmp_b/1e9); %we want to average repetitions at given b-value
    sigs(c_g, 1:numel(b_groups)) = tmp_s; %dont normalise here
end

%find waveforms with fixed gamma and fixed Vomega (code block below is not intended to be read)
vom = gam_vom(:,2);
x_ind = zeros(round(numel(vom)/2), 1);
%look for any waveform on the x-axis
tol = 2; %vomegas this close to each other are considered equal
for c_vom = 1:numel(vom)
    tmp = abs(vom(c_vom) - vom(1:numel(vom)~=c_vom));
    if any(tmp < tol)
        x_ind(1) = c_vom;
        vom_copy = vom(1:numel(vom)~=c_vom);
        vom_on_x = vom_copy(tmp < tol);
        tmp_x_ind = zeros(size(vom_on_x));
        for c_c_vom = 1:numel(vom_on_x)
            tmp_x_ind(c_c_vom) = find(vom == vom_on_x(c_c_vom));
        end
        x_ind(2:end) = tmp_x_ind;
    break
    end
end
%okay, now get y-axis waveforms
y_ind = setdiff((1:numel(vom)), x_ind);
%--------------------------------------------------------------------------


%plot signal vs gamma
sigs = log(sigs);

sp1 = axes(f);%subplot(1,2,1);
for c_b = 1:numel(b_groups)
    x = gam_vom(x_ind, 1); %expecting that these are the waveforms with different gammas
    y = sigs(x_ind, c_b);
    [x, sort_idx] = sort(x, 'ascend');
    y = y(sort_idx);
    plt = plot(sp1, x,y, 'o-', 'MarkerSize', 8, 'LineWidth', 1.5);
    set(plt, 'MarkerFaceColor', plt.Color)
    hold(sp1, 'on')
end
xlabel(sp1, '\Gamma [ms]')
ylabel(sp1, {'\leftarrow Increasing b-value', 'Log Signal'})
grid(sp1, 'minor')
set(sp1, 'MinorGridLineStyle', '-')
% sp1.Position = [0.1500    0.6200    0.7500    0.3000];
%left bottom width height
sp1.Position = [0.3800    0.3200    0.200    0.5000];


%plot signal vs vomega

sp2 = axes(f);%subplot(1,2,2);
for c_b = 1:numel(b_groups)
    x = gam_vom(y_ind, 2); %expecting that these are the waveforms with different vomegas
    y = sigs(y_ind, c_b);
    [x, sort_idx] = sort(x, 'ascend');
    y = y(sort_idx);
    plt2 = plot(sp2, x,y, 'o-', 'MarkerSize', 8, 'LineWidth', 1.5);
    set(plt2, 'MarkerFaceColor', plt2.Color)
    hold(sp2, 'on')
end
xlabel(sp2, ' V_{\omega}^{-1/2} [ms]')
ylabel(sp2, {'\leftarrow Increasing b-value', 'Log Signal'})
grid(sp2, 'minor')
set(sp2, 'MinorGridLineStyle', '-')
sp2.Position =  [0.6500    0.320    0.2000    0.500];



end