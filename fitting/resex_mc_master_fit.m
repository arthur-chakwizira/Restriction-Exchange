%resex_mc_master_fit
%This is the master script which shows how to fit models to simulated data
close all
clc


%SETUP options-----------------------
b_s = (0.001:0.5:2.1)*1e9; %desired b_values
opt = resex_mc_exp_opt(); %xperimental parameters (max gradient amplitude, slew rate, etc)

gwf_fn = 'sde_gwf_for_v4'; %example protocol
%prepare protocol for signal generation
load(gwf_fn);
ngwf = size(gwf, 1);
gwf_fn = 'tmp_gwf';
save(gwf_fn, 'gwf', 'dt')
gwf_fn = 'tmp_gwf'; 

%specify trajectory file path
r_fn = "d_5_k_20_new"; 
%--------------------------------
                                                                               

%FEXI xps (if needed)___________________________
% gwf_fn = "fexi_protocol.mat";
% bs_f = b_s;
% bs_d = b_s;
% xps = resex_mc_xps_fexi_build(gwf_fn, bs_f, bs_d, opt);
%--------------------------------


%FWF xps-----------------------
xps = resex_mc_xps_build(gwf_fn, b_s, opt, r_fn); 
%--------------------------------


%RUN simulation----------------------- (not needed if positions are loaded from file)
% resex_mc_r_save(r_fn, opt)
%--------------------------------


%UNCOMMENT to invoke load positions from different path than given in r_fn
% [file, path] = uigetfile('*.mat','Select file with saved positions.'); %invoke file select dialog
% if file ~= 0 %if user does not cancel dialog
%     r_fn = fullfile(path, file);
% else
%     return
% end 
%--------------------------------


%GENERATE signal using the saved positions and experimental struct
[s_simulated, ~] = resex_mc_signal_from_r(r_fn, xps); %this calculates signal using those saved positions
%--------------------------------


%PREDICT signal using diffnt methods (if you want)
% m =  2e-9;
% opt.exp_opt.using_sim_data = false;
% xps = resex_mc_xps_build(gwf_fn, b_s, opt, r_fn); 
% s_predicted = resex_mc_signal_monoexp(m, xps);
% myfig;
% plot(log(s_simulated), 'go-', 'MarkerFaceColor', 'g');
% plot(log(s_predicted), 'k-')
% legend(["Simulation", "Model"])


%FIT to simulated data (uncomment correct method)-----------
% mfs = resex_mc_fit_fm_any_geo(s_simulated, xps) %full model, any geometry
% mfs = resex_mc_fit_fexi(s_simulated, xps); %fexi
% mfs = resex_mc_fit_myrep(s_simulated, xps) %cumulant expasions (direct cumulants)
% mfs = resex_mc_fit_myrep_new(s_simulated, xps) %cumulant expasions (direct cumulants)
% mfs = resex_mc_fit_fm(s_simulated, xps) %full analytical 2-compartment model
% mfs = resex_mc_fit_fm_new(s_simulated, xps); %this is just for fun, AND IT DOESN'T WORK
mfs = resex_mc_fit_rep(s_simulated, xps) %representation from ABSTRACT, Nilsson 2019
% mfs = resex_mc_fit_rep_simple(s_simulated, xps) %representation from ABSTRACT but with V_D, C_DR, V_R bunched up into V
% mfs = resex_mc_fit_monoexp(s_simulated, xps) %monoexponential fit
% mfs = resex_mc_fit_km(s_simulated, xps) %modified Kärger model
%--------------------------------


%DISPLAY goodness of fit-----------
figure('Color', 'w')
plot(s_simulated, 'bo', 'MarkerFaceColor', 'b')
hold on
plot(mfs.signal, 'rs-', 'LineWidth', 1, 'MarkerEdgeColor', 'r')
legend({'Simulated', 'Model Fit'})
set(gca,'tickdir', 'out', 'ticklength', [0.01 0.01], 'linewidth', 1.5, 'fontsize', 20, 'fontname', 'arial', 'box', 'on');
ylabel('Signal')
title('GOODNESS OF FIT')
%--------------------------------



%PRECISION study: uncomment to add noise and fit to all samples in loop
% snr0 = 200;
% T2 = 80;
% snr = snr0;%*exp(-TE/T2);
% N_samp = 100;
% s_noise = aux_add_noise(s_simulated, snr, N_samp);
% 
% pars = zeros(N_samp, 6); %f1, f2, D1, D2, d, k
% sigs = zeros(N_samp, numel(s_simulated));
% for c_s = 1:N_samp
%    s_sim = s_noise(c_s, :);
%    tmp_pars = resex_mc_fit_rep(s_sim', xps);
%    pars(c_s, :)  = tmp_pars.params;
%    sigs(c_s, :) = tmp_pars.signal;
%    clc
%    disp(c_s)
% end
% hold off
% histogram(pars(:,end))