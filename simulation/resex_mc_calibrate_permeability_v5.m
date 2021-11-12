%calibrate exchange-permeability automatically
%its time-consuming to do this manually
%We want a general relationship between input permeability and exchange
%rate
%we will use a dt of 1e-5 in all cases

%%

r_fn = "junk.mat"; %path to save positions to
%--------------------------------                                                                    



ds = [2 4 5 6 8 10 12 14 16 18 20 25 30 40 50 60 70 80 90 100]*1e-6;
kappas = [1e-6  2e-6 5e-6 10e-6 20e-6 100e-6 200e-6 300e-6 400e-6 500e-6];

D_1 = 2e-9;
D_2 = 2e-9;
calib.D_1 = D_1;
calib.D_2 = D_2;


%% calibrate k_12
fprintf(1,'Running simulation : %3d%%\n',0);
prog = 0;

kappa_21 = 0;
f_1 = 1;
for c_d = 1:numel(ds)

    d = ds(c_d);
    calib.kappa_12.("d_" + num2str(d*1e6)) = zeros(numel(kappas), 2);
    for c_kap = 1:numel(kappas)
        kappa_12 = kappas(c_kap);
    opt = calibration_sim_opt(opt, kappa_12, kappa_21, d, D_1, D_2, f_1); %Simulation options 
    resex_mc_r_save(r_fn, opt)
    k = calibrate_get_k_from_n(1);
    calib.kappa_12.("d_" + num2str(d*1e6))(c_kap, :) = [kappas(c_kap)  k];
    
        prog = prog+1;
        fprintf(1,'\b\b\b\b%3.0f%%', 100*(prog/(2*numel(ds)*numel(kappas))) );
    end
end

%% calibrate k_21
kappa_12 = 0;
f_1 = 0;
for c_d = 1:numel(ds)
   
    d = ds(c_d);
    calib.kappa_21.("d_" + num2str(d*1e6)) = zeros(numel(kappas), 2);
    for c_kap = 1:numel(kappas)
        kappa_21 = kappas(c_kap);
    opt = []; %simulation options
    opt = calibration_sim_opt(opt, kappa_12, kappa_21, d, D_1, D_2, f_1); %Simulation options 
    resex_mc_r_save(r_fn, opt)
    k = calibrate_get_k_from_n(2);
    calib.kappa_21.("d_" + num2str(d*1e6))(c_kap, :) = [kappas(c_kap)  k];
    
            prog = prog+1;
        fprintf(1,'\b\b\b\b%3.0f%%', 100*(prog/(2*numel(ds)*numel(kappas))) );
    end
end
%--------------------------------
 save('resex_mc_calib_v5', 'calib')
%%
function k = calibrate_get_k_from_n(popn_id)
if popn_id == 1
clear 'popn_1' 'time'
load('resex_mc_popn_1.mat', 'popn_1', 'time')
N = popn_1(popn_1~=0);
t = time(popn_1~=0);
else
  clear 'popn_2' 'time'
load('resex_mc_popn_1.mat', 'popn_2', 'time')
N = popn_2(popn_2~=0);
t = time(popn_2~=0);  
end
% N = popn_1(1:1000);
% t = time(1:1000);
logN = log(N);
p = polyfit(t', logN, 1);   

k = abs(p(1));

% figure
% set(gcf, 'color', 'w');
% set(gca,'box','on','linewidth',0.5,'layer','top', 'FontWeight', 'normal', 'FontSize', 12)
% xlabel('time [s]')
% ylabel('Intracellular population (log)')
% hold on
% tf = linspace(min(t), max(t), 1000);
% Nf = polyval(p,tf);
% plot(tf, Nf, 'LineWidth', 2)
% plot(t,logN,'--', 'LineWidth',2)
% legend({'Fit', 'Data'})

end


function opt = calibration_sim_opt(opt, kappa_12, kappa_21, d, D_1, D_2, f_1)
if isempty(opt); opt.present = 1; end
opt.sim_opt.D_1 = D_1; %intracellular diffusivity
opt.sim_opt.D_2 = D_2; %extracellular diffusivity


opt.sim_opt.f_1 = f_1;
opt.sim_opt.f_2 = 1-opt.sim_opt.f_1;


opt.sim_opt.kappa_12 = kappa_12; %permeability
opt.sim_opt.kappa_21 = kappa_21; %permeability

opt.sim_opt.r = (d/2); %intracellular radius



opt.sim_opt.N_part = 1e3; %total population


opt.sim_opt.T = 100e-3; %total simulation time
opt.sim_opt.sim_dt = 1e-4; %simulation resolution
opt.sim_opt.do_avg = false; %average positions
opt.sim_opt.do_samp = true; %average positions
opt.sim_opt.avg_dt = 1e-3; %average every millisecond;

v_2 = sqrt(4*opt.sim_opt.D_2/opt.sim_opt.sim_dt); %particle velocity in 2. This is adapted from Szafer A, Zhong J, Gore JC. 1995 and also Meier C, et al. 2003
opt.sim_opt.p_21 = 4*opt.sim_opt.kappa_21/v_2; %kappa here is the membrane permeability
v_1 = sqrt(4*opt.sim_opt.D_1/opt.sim_opt.sim_dt); %particle velocity in 1. This is adapted from Szafer A, Zhong J, Gore JC. 1995 and also Meier C, et al. 2003
opt.sim_opt.p_12 = 4*opt.sim_opt.kappa_12/v_1; %kappa here is the membrane permeability

opt.sim_opt.n_dim = 2; %doing simulations in 2D plane

% %get required number of cells
% opt.sim_opt.N_cells = resex_mc_get_num_cells(opt.sim_opt.sim_dt, opt.sim_opt.T, ...
%     max(opt.sim_opt.D_1, opt.sim_opt.D_2), 2*opt.sim_opt.r, opt.sim_opt.n_dim);
% if mod(ceil(opt.sim_opt.N_cells/2), 2) == 0; errordlg('Selected number of cells is invalid'); end

opt.sim_opt.loiter = 100;

%visualisation of simulation
opt.sim_opt.show = false; %display initial configuration or not

opt.sim_opt.film.make = false; %disables filming if false
opt.sim_opt.film.save_film = false;
opt.sim_opt.film.film_title = 'resex_mc_film.avi';
opt.sim_opt.film.frame_rate = 20;

end
