
%This is the master script which controls the entire simulation framework
close all



%SETUP options-----------------------
opt = resex_mc_exp_opt(); %experimental options
% 
% ds = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20];
% ds = flip(ds);
% ks = [0 2 4 6 8 10 12 14 16 18 20];
% ks = flip(ks);

%start with most important case
ds = 10;
ks = 10;
for c_d = 1:numel(ds)
    opt.tmp.d = ds(c_d)*1e-6;
    for c_k = 1:numel(ks)
        opt.tmp.k = ks(c_k);
        opt = resex_mc_sim_opt(opt, "two_exchanging_compartments"); %Simulation options
% 
% r_fn = 'junk.mat';
r_fn = "d_" + num2str(ds(c_d)) + "_k_" + num2str(ks(c_k)) + "_v3.mat"; %v3 means that these were generated using ReSex_MC_v3 functions
%--------------------------------

if ~isfile(r_fn) %don't want to repeat
%RUN simulation-----------------------
resex_mc_r_save(r_fn, opt)
%--------------------------------
end
    end
end

