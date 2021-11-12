%waveform_generator_for_gamma_vomega_dde
%Generates DDE waveforms

clc
close all
tic
%% Define waveform parameters
% wps.T = 100e-3; %total encoding time in seconds
wps.dt = 1e-4; %time-step in seconds
wps.g_max = 75e-3; %maximum gradient amplitude in Tesla
wps.slew_max = 71; %maximum slew rate in Tesla/seconds
wps.slew_rate = 70;
wps.t_start = 0;
wps.gam_ma = 2*pi*42.6e6; %gamma constant for the 1H nucleus, 42.6e6 Hz/T

g_max = wps.g_max;

%%
Nwfs = 10; %number of waveforms
Ntp = 3001; %number of time points
wfs = zeros(Nwfs, Ntp);

qs = zeros(size(wfs));
b_values = zeros(Nwfs,1);
Ns = zeros(Nwfs,1); %number of steps in waveform
gammas = zeros(Nwfs, 1);
vomegas = zeros(Nwfs, 1);
TEs = zeros(Nwfs, 1);

gam_ma = wps.gam_ma;
dt = wps.dt;

%%
wf_tracker = zeros(1, size(wfs,2));
trash = 0;

b = 5e9;
bmax = 10e9;
min_sep = 4e-3;

fix_b = true;

ef = zeros(Nwfs, 1);


for c_wf = 1:Nwfs
    reject_wf = true;
    while reject_wf
        
        wps.T = 20e-3+ randi([0 230])*1e-3; %total encoding time in seconds
        t_max= (floor(rand()*wps.T*(1/1e-3)/2)-4)*1e-3; %49e-3; %waveform ends here
        wps.time = linspace(0, t_max, round(t_max/wps.dt)+1); %time vector = linspace(0, T, N_steps) in seconds
        wps.N_steps = length(wps.time); %number of time-steps = round(T/dt)+1
        
        try
        wps.sde.delta = randi([5 floor((t_max-min_sep)*0.5*(1e3))])*1e-3; %gradient pulse +duration in seconds
        catch
            continue
        end
        
        wps.sde.interval = t_max - 2*wps.sde.delta-1e-3;
                
        wps.g_max = rand()*g_max; %allow freedom of height
        
        
        wf_info = waveform_sde(wps);
        
        tmp_time  = linspace(0, wps.T, round(wps.T/wps.dt));
        tmp_g = wf_info.waveform';
        g = zeros(1, numel(tmp_time));
        g(1:numel(tmp_g)) = tmp_g;
        g((end-numel(tmp_g)+1):end) = tmp_g;
        
        q = msf_const_gamma()*cumsum(g)*dt;
        b_value =sum(q.^2)*dt;
                
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
        
        
        TE = numel(g)*dt -dt;
                
 if TE*1e3 > 250; reject_wf = true; end
        
    end
        eff2 = 4*b_value/(wps.gam_ma^2*max(g)^2*wps.T^3); %efficiency

    %%
    
    
    wfs(c_wf,1:numel(g)) = g;
    qs(c_wf, 1:numel(q)) = q;
    b_values(c_wf) = b_value;
    Ns(c_wf) = wps.N_steps;
    gammas(c_wf) = 1e3*gam;
    vomegas(c_wf) = 1e3/sqrt(vom);
    TEs(c_wf) = TE*1e3;
        ef(c_wf) = eff2;
    clc
    disp(strcat('Generated :', num2str(c_wf), ' of :', num2str(Nwfs)))
end

% save('/Users/arthur/Documents/LUND_UNIVERSITY/PHD/SIMULATION_CODE/cpi/dde_landscape_b5.mat', 'wfs','qs','Ns','b_values','bmax', 'wps', 'gammas', 'vomegas', 'TEs', 'ef')


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
