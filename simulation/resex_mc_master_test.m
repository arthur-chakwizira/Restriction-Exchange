

%This is the master script which controls the entire simulation framework
%V5 is similar to V4 except it:
%   reverts to fixed jump lengths in random directions
%   fixes the spacing between cells in substrate
%   
%V4 of the code attempts to improve performance by rearraging data to
%ensure looping column by column as much as possible. This means that
%positions are saved as N_particles x N_dimensions x N_time_points which is
%not compatible with previous versions of this code.
%You can call resex_mc_reshape_trajectory to reshape old trajectory files
%to work with this new code.
%In V4, time is almost always the last dimension.
%V3 is compatible with all past versions of the code, with the difference
%that it uses boundary conditions for the diffusion process
close all
clc

%SETUP options-----------------------
opt = resex_mc_exp_opt(); %experimental options

%--------------------------------


% RUN simulation-----------------------
key = "two_exchanging_pools";
base_folder = "/Users/arthur/Documents/LUND_UNIVERSITY/PHD/PAPER_3/v4/TRAJECTORIES/";
r_fn = base_folder + key;

% % opt.tmp.d = 10e-6;
% % opt.tmp.k = 0;
% % opt = resex_mc_sim_opt(opt, "two_exchanging_compartments"); %Simulation options
% % resex_mc_r_save(r_fn, opt)

% r_fn = 'one_gaussian_pool';
% opt = resex_mc_sim_opt(opt, "one_gaussian_pool"); 
% resex_mc_one_gaussian_pool_r_save(r_fn, opt)

% r_fn = 'three_gaussian_pools';
% opt = resex_mc_sim_opt(opt, "three_gaussian_pools"); 
% resex_mc_three_gaussian_pools_r_save(r_fn, opt)
% 

% r_fn = 'one_circle';
% opt = resex_mc_sim_opt(opt, "one_circle"); 
% resex_mc_one_circle_r_save(r_fn, opt)

% r_fn = 'two_non_exchanging_pools'; %'two_exchanging_pools';
% opt = resex_mc_sim_opt(opt, "two_exchanging_pools"); 
% resex_mc_two_exchanging_pools_r_save(r_fn, opt)

% r_fn = 'two_exchanging_pools'; %'two_exchanging_pools';
opt = resex_mc_sim_opt(opt, "two_exchanging_pools"); 
resex_mc_two_exchanging_pools_r_save(r_fn, opt)
