%Random generation of waveforms of type sine-modulated OGSE
tic
close all
% dock_figures(1)
%% Define waveform parameters
wps.dt = 1e-4; %time-step in seconds
wps.g_max = 75e-3; %maximum gradient amplitude in Tesla
wps.slew_max = 70; %maximum slew rate in Tesla/seconds
wps.interval = 1e-4;%DELTAs(c_d);%intervals(c_D); %time interval between end of pulse 1 and start of pulse 2 in seconds (interval + delta = DELTA)
wps.t_start = 0;%floor(0.5*(wps.T - ...
wps.gam_ma = 2*pi*42.6e6; %gamma constant for the 1H nucleus, 42.6e6 Hz/T

%% generate waveforms
Nwfs = 100; %number of waveforms
Ntp = 3001; %number of time points
wfs = zeros(Nwfs, Ntp);
qs = zeros(size(wfs));
b_values =  zeros(Nwfs,1);
gammas = zeros(Nwfs, 1);
vomegas = zeros(Nwfs, 1);
Ns = zeros(Nwfs,1); %number of steps in waveform
gam_ma = wps.gam_ma;

b = 5e9; %generate waveforms at this b value
bmax = 10e9;
min_sep = 4e-3; %ms

g_max = wps.g_max;
dt = wps.dt;

TEs = zeros(Nwfs, 1);
ef = zeros(Nwfs, 1);

fix_b = 0; %demand fixed b-value or not

for c_wf =  1:Nwfs
    reject_wf = true;
    reject_cp = true;
    while reject_wf
        wps.T = 15e-3+ randi([0 235])*1e-3; %total encoding time in seconds %
        t_max= (floor((wps.T*1e3)/2)-1)*1e-3; %pulse ends here
        t_pulse = 0:dt:t_max;
        wps.delta = t_max; %gradient pulse +duration in seconds
        
        t_end = rand*(t_max-min_sep);
        
        t_lobe = 0:dt:t_end;
        
        n_lobes = randi([2 40]); %quicker than above
        
        
        %get pulse
        pulse = t_pulse*0;
        period = 2*max(t_lobe)/n_lobes;
        tmp_pulse = wps.g_max*sin((2*pi/period)*t_lobe);
        pulse(1:numel(tmp_pulse)) = tmp_pulse;
        %force maximum amplitude
        %         pulse = wps.g_max*pulse./max(pulse);
        
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
        
        %impose some distanc between points
        %                 mn = [(gammas(~~gammas)), (vomegas(~~vomegas))];
        %                 d = (1e3*gam - mn(:, 1)).^2 + (1e3/sqrt(vom) - mn(:, 2)).^2;
        %                 if ~isnan(d) & any(d<1); reject_wf = true; end
        
        
        if ismember([gam*1e3  1e3/sqrt(vom)], [gammas(1:c_wf) vomegas(1:c_wf)], 'rows')
            reject_wf = true; end
        
        TE = numel(g)*dt -dt;        
        if TE*1e3 > 250; reject_wf = true; end
        
    end
    
    wfs(c_wf,1:wps.N_steps) = g;
    qs(c_wf, 1:wps.N_steps) = q;
    b_values(c_wf) = b_value;
    Ns(c_wf) = wps.N_steps;
    gammas(c_wf) = 1e3*gam;
    vomegas(c_wf) = 1e3/sqrt(vom);
    TEs(c_wf) = TE*1e3;
    
    clc
    disp(strcat('Generated :', num2str(c_wf), ' of :', num2str(Nwfs)))
end




save('ogse_sine_gwf_example.mat', 'wfs','qs','Ns','b_values','bmax', 'wps', 'gammas', 'vomegas', 'TEs')


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
numj = min([numel(j), 4]);
numcols = 4;
numrows = ceil(numj/numcols);
for c_j = 1:numj
    subplot(numrows, numcols, c_j)
    fill(1:size(wfs,2), wfs(j(c_j),:)*1e3, 'b')
    title(num2str(j(c_j)))
    ylabel('g [mT/m]')
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
