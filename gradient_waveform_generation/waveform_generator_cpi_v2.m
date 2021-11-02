%
%FWF waveform generator

tic
close all
dock_figures(1)
%% Define waveform parameters
% wps.T = 100e-3; %total encoding time in seconds
wps.dt = 1e-4; %time-step in seconds
wps.g_max = 75e-3; %maximum gradient amplitude in Tesla
wps.slew_max = 70; %maximum slew rate in Tesla/seconds
% wps.delta = 49e-3; %gradient pulse duration in seconds
wps.interval = 1e-4;%DELTAs(c_d);%intervals(c_D); %time interval between end of pulse 1 and start of pulse 2 in seconds (interval + delta = DELTA)
wps.t_start = 0;%floor(0.5*(wps.T - ...x
% 2*wps.delta-wps.interval)/wps.dt)*wps.dt;
% wps.time = linspace(0, wps.T,round(wps.T/wps.dt)+1); %time vector = linspace(0, T, N_steps) in seconds
% wps.N_steps = length(wps.time); %number of time-steps = round(T/dt)+1
wps.gam_ma = 2*pi*42.6e6; %gamma constant for the 1H nucleus, 42.6e6 Hz/T

%% generate waveforms
Ts = (90:250)*1e-3;
n_eff_checks = 10;
n_g_per_te = 100;


Nwfs = n_g_per_te*numel(Ts); %number of waveforms
Ntp = 3001; %number of time points
wfs = zeros(Nwfs, Ntp);
qs = zeros(size(wfs));
b_values =  zeros(Nwfs,1);
gammas = zeros(Nwfs, 1);
vomegas = zeros(Nwfs, 1);
Ns = zeros(Nwfs,1); %number of steps in waveform
gam_ma = wps.gam_ma;
% dt = wps.dt;
% cps = zeros(Nwfs, 5); %control points

wf_tracker = zeros(1, size(wfs,2));
trash = 0;
b = 5e9; %generate waveforms at this b value
bmax = 10e9;
min_sep = 4e-3; %ms

g_max = wps.g_max;
dt = wps.dt;

TEs = zeros(Nwfs, 1);
ef = zeros(Nwfs, 1);


fix_b = 1;
n_points = 7;
TE_thres = 5;


c = 0;

for c_T =  1:numel(Ts)
    wps.T = Ts(c_T);
    for c_wf = 1:n_g_per_te
        eff = 0; %choose most efficient of 100 trials
        %temporary arrays for choosing most efficient with longest Y
        tmp_eff = zeros(n_eff_checks, 1);
        tmp_g = zeros(n_eff_checks, Ntp);
        tmp_q = zeros(n_eff_checks, Ntp);
        tmp_b = zeros(n_eff_checks, 1);
        tmp_N_steps =  zeros(n_eff_checks, 1);
        tmp_gam = zeros(n_eff_checks, 1);
        tmp_vom = zeros(n_eff_checks, 1);
        tmp_TE = zeros(n_eff_checks, 1);
        
%         tmp_power_sigma = zeros(n_eff_checks, 1);
            
        for trial = 1:n_eff_checks
            reject_wf = true;
            while reject_wf
                t_max= (floor((wps.T*1e3)/2)-1)*1e-3; %pulse ends here
                t_pulse = 0:dt:t_max;
                wps.delta = t_max; %gradient pulse +duration in seconds
                
                t_end = rand*(t_max-min_sep);
                
                
                
                t_lobe = 0:dt:t_end;
                
                %get cp
                [x,y] = get_cpi(t_end, g_max);
                if x == 0; reject_wf = true; continue; end
                
                tmp_pulse = interp1(x,y, t_lobe, 'pchip');
                tmp_pulse(abs(tmp_pulse) < 0.5e-3) = 0;
                
                
                
                pulse = t_pulse*0;
                pulse(1:numel(tmp_pulse)) = tmp_pulse;
                
                
                %force maximum amplitude
                pulse = wps.g_max*pulse./max(abs(pulse));
                
                %%
                wps.pulse = pulse;
                
                wps.time = linspace(0, wps.T,round(wps.T/wps.dt)+1); %time vector = linspace(0, T, N_steps) in seconds
                wps.N_steps = length(wps.time); %number of time-steps = round(T/dt)+1
                
                wf_info = get_wf_from_pulse(wps);
                g = wf_info.waveform;
                b_value = wf_info.b_value;
                q = wf_info.q_vector;
                
                
                %g is not the same length as time vector anymore due to cropping
                wps.N_steps = numel(g);
                
                
                if fix_b
                    g_scale = sqrt(b/b_value);
                    g = g_scale*g;
                    q = gam_ma*cumsum(g)*dt;
                    b_value = sum(q.^2)*dt;
                end
                
                reject_wf = reject_pulse(g, wps);
                
                [gam, vom] = resex_mc_gwf_to_gamma_vomega(g,q,b_value,dt);
                
                
                %                 mn = [(gammas(~~gammas)), (vomegas(~~vomegas))];
%                 d = (1e3*gam - mn(:, 1)).^2 + (1e3/sqrt(vom) - mn(:, 2)).^2;
%                 if ~isnan(d) & any(d<1); reject_wf = true; end
                
                
                
                %         if ismember([gam*1e3  1e3/sqrt(vom)], [gammas(1:c_wf) vomegas(1:c_wf)], 'rows')
                %             reject_wf = true; end
                
                %reject waveforms whose TE/Gamma is worse than best SDE case (3)
                [N, g] = N_from_gwf(g);
                TE = (N-1)*dt;

                %         if TE > TE_thres*gam; reject_wf = true; end
                         if TE*1e3 > 250; reject_wf = true; end

                
            end
            
             eff2 = 4*b_value/(wps.gam_ma^2*max(g)^2*wps.T^3); %efficiency
%              [Qw, w] = resex_mc_qt_to_qw(q, dt);
%              fm = sum(Qw(w>0).*w(w>0))*(w(end)-w(end-1));
%              power_std = fm;
             
        tmp_eff(trial) = eff2;
        tmp_g(trial, 1:numel(g)) = g;
        tmp_q(trial, 1:numel(q)) = q;
        tmp_b(trial) = b_value;
        tmp_N_steps(trial) = wps.N_steps;
        tmp_gam(trial) = gam;
        tmp_vom(trial) = vom;
        tmp_TE(trial) = TE;
        
%         tmp_power_sigma(trial) = power_std;
           
        end
        %now we will choose the best of those waveforms
        c = c+1;
        
        metric = tmp_eff.*tmp_gam.*tmp_vom;
%         metric = tmp_power_sigma;
        
        [~, ind] = max(metric);
        
        
        wfs(c,:) = tmp_g(ind, :);
        qs(c, :) = tmp_q(ind, :);
        b_values(c) = tmp_b(ind);
        Ns(c) = tmp_N_steps(ind);
        gammas(c) = 1e3*tmp_gam(ind);
        vomegas(c) = 1e3/sqrt(tmp_vom(ind));
        TEs(c) = tmp_TE(ind)*1e3;
        ef(c) = tmp_eff(ind);
        clc
        disp(strcat('Generated :', num2str(c), ' of :', num2str(Nwfs)))
    end
end




% save('/Users/arthur/Documents/LUND_UNIVERSITY/PHD/SIMULATION_CODE/cpi/try_efficiency_constraint_new.mat', 'wfs','qs','Ns','b_values','bmax', 'wps', 'gammas', 'vomegas', 'TEs', 'ef')
% save('/Users/arthur/Documents/LUND_UNIVERSITY/PHD/SIMULATION_CODE/cpi/fwf_landscape_b5.mat', 'wfs','qs','Ns','b_values','bmax', 'wps', 'gammas', 'vomegas', 'TEs', 'ef')


%
figure; hold on
plot(gammas, vomegas, 'o', 'MarkerFaceColor', [0.5 0.5 0.5])
j = boundary(gammas,vomegas,0.1);
hold on;
plot(gammas(j),vomegas(j), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlabel('$\Gamma$ [ms]', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$1/\sqrt V_\omega$ [ms]',  'Interpreter', 'latex', 'FontSize', 30)
set(gca, 'box', 'on', 'Color', 'w', 'FontSize', 20)


figure
numj = numel(j);
numcols = 4;
numrows = min([4, ceil(numj/numcols)]);
for c_j = 1:numj
    subplot(numrows, numcols, c_j)
    fill(1:size(wfs,2), wfs(j(c_j),:), 'b')
    title(num2str(j(c_j)))
end

%plot T vs Gamma
figure()
plot(gammas, TEs, 'ro', 'MarkerFaceColor', [1 0.4 0.4], 'HandleVisibility', 'off')
hold on
plot(linspace(min(gammas), max(gammas)), 3*linspace(min(gammas), max(gammas)), ...
    'k-', 'LineWidth', 2)
xlabel('$\Gamma$ [ms]','Interpreter', 'latex', 'FontSize', 30)
ylabel('$T$ [ms]', 'Interpreter', 'latex', 'FontSize', 30)
set(gca, 'box', 'on', 'Color', 'w', 'FontSize', 20)
legend("$T = 3 \cdot \Gamma$",'Interpreter', 'latex', 'FontSize', 20)
toc


function [x,y] = get_cpi(t_end, g_max)
%generate cpi
min_dt = 2e-3; %we want at least 3ms between consecutive cpi

n = floor(t_end/min_dt); %how many points, max
x = NaN(n, 1);
y = NaN(n, 1);
x(1) = 0; y(1) = 0;
for c_n = 2:n
    if any(x>t_end); break; end
    x(c_n) = x(c_n-1) + min_dt + (t_end-x(c_n-1) - min_dt)*rand();
    y(c_n) = -g_max+ 2*rand()*g_max;
end
x(isnan(x)) = [];
y(isnan(y)) = [];
y(end) = 0;
if numel(y) > 2; y(end-1) = y(end-2)*rand; end
end
