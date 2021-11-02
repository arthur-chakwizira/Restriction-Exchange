%waveform_generator_for_gamma_vomega_sde
%adding efficiency constraint
%Also trying to limit power to low frequencies

tic
close all
% dock_figures(1)
%% Define waveform parameters
wps.dt = 1e-4; %time-step in seconds
wps.g_max = 75e-3; %maximum gradient amplitude in Tesla
wps.slew_max = 70; %maximum slew rate in Tesla/seconds
wps.interval = 1e-4;%DELTAs(c_d);%intervals(c_D); %time interval between end of pulse 1 and start of pulse 2 in seconds (interval + delta = DELTA)
wps.t_start = 0;%floor(0.5*(wps.T - ...x
wps.gam_ma = 2*pi*42.6e6; %gamma constant for the 1H nucleus, 42.6e6 Hz/T

%% generate waveforms
Ts = (50:250)*1e-3;
n_eff_checks = 20;

Nwfs = 20100;%n_g_per_te*numel(Ts); %number of waveforms
Ntp = 3001; %number of time points
wfs = zeros(Nwfs, Ntp);
qs = zeros(size(wfs));
b_values =  zeros(Nwfs,1);
gammas = zeros(Nwfs, 1);
vomegas = zeros(Nwfs, 1);
Ns = zeros(Nwfs,1); %number of steps in waveform
gam_ma = wps.gam_ma;

wf_tracker = zeros(1, size(wfs,2));
trash = 0;
b = 2e9; %generate waveforms at this b value
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
    if wps.T <= 50e-3; n_g_per_te = 9; end
    if wps.T > 50e-3 && wps.T < 60e-3; n_g_per_te = 10; end
    if wps.T >= 60e-3 && wps.T < 100e-3; n_g_per_te = 30; end
    if wps.T >= 100e-3; n_g_per_te = 100; end
    
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
                    
        for trial = 1:n_eff_checks
            reject_wf = true;
            while reject_wf
                t_max= (floor((wps.T*1e3)/2)-1)*1e-3; %pulse ends here
                wps.time = linspace(0, wps.T,round(wps.T/wps.dt)+1); %time vector = linspace(0, T, N_steps) in seconds
                wps.N_steps = length(wps.time); %number of time-steps = round(T/dt)+1
                
                wps.sde.delta = randi([5 floor((t_max-min_sep)*(1e3))])*1e-3; %gradient pulse +duration in seconds
                
                wps.sde.interval = wps.T - 2*wps.sde.delta-1e-3;
                
                wps.g_max = rand()*g_max; %allow freedom of height
                
                wf_info = waveform_sde(wps);
                
               g = wf_info.waveform';
               b_value = wf_info.b_value;
               q = wf_info.q_vector';
 
                
                if fix_b
                    g_scale = sqrt(b/b_value);
                    g = g_scale*g;
                    q = gam_ma*cumsum(g)*dt;
                    b_value = sum(q.^2)*dt;
                end
                
                reject_wf = reject_pulse(g, wps);
                
                [gam, vom] = resex_mc_gwf_to_gamma_vomega(g,q,b_value,dt);
                
                
                        if ismember([gam*1e3  1e3/sqrt(vom)], [gammas(1:c_wf) vomegas(1:c_wf)], 'rows')
                            reject_wf = true; end
                
                %reject waveforms whose TE/Gamma is worse than best SDE case (3)
                TE = numel(g)*dt -dt;
                
            end
            
             eff2 = 4*b_value/(wps.gam_ma^2*max(g)^2*wps.T^3); %efficiency
             
        tmp_eff(trial) = eff2;
        tmp_g(trial, 1:numel(g)) = g;
        tmp_q(trial, 1:numel(q)) = q;
        tmp_b(trial) = b_value;
        tmp_N_steps(trial) = wps.N_steps;
        tmp_gam(trial) = gam;
        tmp_vom(trial) = vom;
        tmp_TE(trial) = TE;
                   
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


garbage = (gammas == 0 & vomegas == 0);
wfs(garbage, :) = [];
qs(garbage, :) = [];
Ns(garbage) = [];
b_values(garbage) = [];
gammas(garbage) = [];
vomegas(garbage) = [];
TEs(garbage) = [];
ef(garbage) = [];


% % save('/Users/arthur/Documents/LUND_UNIVERSITY/PHD/SIMULATION_CODE/cpi/try_optimise_sde.mat', ...
%     'wfs','qs','Ns','b_values','bmax', 'wps', 'gammas', 'vomegas', 'TEs', 'ef')


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

%plot TE vs Gamma
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