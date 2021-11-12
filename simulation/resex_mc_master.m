
%This is the master script which controls the simulation framework
%Restriction-exchange simulations on a 2D plane.

%SETUP options-----------------------
opt.tmp.d = 10e-6; %size
opt.tmp.k = 20; %exchange rate (NOTE! you only get this k when the simulation
%               time step and D0 match the settings in the calibration file
%               Consult resex_mc_calibrate_permeability_v5)
%               Exchange rate can also be adjusted by manually setting the
%               permeability (kappa_12, kappa_21) in resex_mc_sim_opt
opt = resex_mc_sim_opt(opt, "two_exchanging_compartments"); %Simulation options
r_fn = "junk"; %"example_simulation" ; %trajectory file path


% RUN simulation-----------------------
resex_mc_r_save(r_fn, opt)