
%This is the master script which controls the entire simulation framework
% close all
clc


%SETUP options-----------------------
b_s = (0.001:0.5:5.1)*1e9; %desired b_values
% parent_path = "/Users/arthur/Documents/LUND_UNIVERSITY/PHD/SIMULATION_CODE/cpi/gamma_vomega_not_enough/28_01_2021_efficiency_consraint/new_fit_results/fitting_without_gam_vom/";
% gwf_fn = parent_path + "protocol_b_2_TE_132.mat"; %path to saved protocol
% 
% gwf_fn = 'TE_120_b5_final_protocol';
gwf_fn = 'gwf_for_vac';
% % gwf_fn = "/Users/arthur/Documents/LUND_UNIVERSITY/PHD/SIMULATION_CODE/cpi/gamma_vomega_not_enough/28_01_2021_efficiency_consraint/new_fit_results/four_waveforms_b_2/protocol_b_2_TE_224";
% %'markus_sde.mat';
load(gwf_fn);
ngwf = 1;
gwf = gwf(1:ngwf, :);
gwf_fn = 'tmp_gwf';
save(gwf_fn, 'gwf', 'dt')
gwf_fn = 'tmp_gwf'; 

opt = [];
% r_fn = 'one_gaussian_pool';
% r_fn = 'two_non_exchanging_comp';
r_fn = "two_non_exchanging_comp";
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


%RUN simulation-----------------------
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


%predict signal using diffnt methods
% m = [1e-9, 2e-9, 0.7, 0];
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
% mfs = resex_mc_fit_rep(s_simulated, xps) %representation from ABSTRACT
% mfs = resex_mc_fit_rep_simple(s_simulated, xps) %representation from ABSTRACT but with V_D, C_DR, V_R bunched up into V
mfs = resex_mc_fit_monoexp(s_simulated, xps) %monoexponential fit
% mfs = resex_mc_fit_km(s_simulated, xps) %modified Kärger model
%--------------------------------

return

%DISPLAY goodness of fit-----------
figure
hold on
plot(s_simulated, 'o')
plot(mfs.signal, 's-')
%--------------------------------

return

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
