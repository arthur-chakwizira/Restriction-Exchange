function opt = resex_mc_sim_opt(opt, env)
%env  = environment
if nargin < 1; opt.present = 1; end

    %general settings
    opt.sim_opt.T = 100e-3; %total simulation time
    opt.sim_opt.sim_dt = 1e-4; %simulation resolution
    opt.sim_opt.do_avg = false; %average positions
    opt.sim_opt.do_samp = true; %sample positions every opt.sim_opt.avg_dt
    opt.sim_opt.avg_dt = 1e-3; %average/sample every opt.sim_opt.avg_dt;
    
    opt.sim_opt.loiter = 500; %delay recording for this many steps
    opt.sim_opt.N_part = 1e2; %total population
    opt.sim_opt.n_dim = 2; %doing simulations in 2D plane
    
    %visualisation of simulation
    opt.sim_opt.show = true; %display initial configuration or not
    
    opt.sim_opt.film.make = false; %disables filming if false
    opt.sim_opt.film.save_film = false;
    opt.sim_opt.film.film_title = 'resex_mc_film_one_particle.avi';
    opt.sim_opt.film.frame_rate = 20;
    
    %specific settings
    if strcmp(env, "two_exchanging_compartments")
        opt.sim_opt.D_1 = 2e-9; %intracellular diffusivity
        opt.sim_opt.D_2 = 2e-9; %extracellular diffusivity
        %these are not in SI units
        d = opt.tmp.d;
        k_in = opt.tmp.k;
        D = opt.sim_opt.D_1;
        
        opt.sim_opt.f_1 = 0.7;
        opt.sim_opt.f_2 = 1-opt.sim_opt.f_1;
        
        k_12 = k_in*(1 - opt.sim_opt.f_1);
        k_21 = k_in*(opt.sim_opt.f_1);
        
        %set permeability manually or quiry calibrated values
        %example: opt.sim_opt.kappa_12 = 3e-4; opt.sim_opt.kappa_21 = 7e-4
        %that gives fast exchange
        opt.sim_opt.kappa_12 = resex_mc_get_kappa_from_k(d, D, k_12, "from_1_to_2"); %permeability
        opt.sim_opt.kappa_21 = resex_mc_get_kappa_from_k(d, D, k_21, "from_2_to_1"); %permeability
        
        opt.sim_opt.r = (d/2); %intracellular radius
        
        v_2 = sqrt(4*opt.sim_opt.D_2/opt.sim_opt.sim_dt); %particle velocity in 2. This is adapted from Szafer A, Zhong J, Gore JC. 1995 and also Meier C, et al. 2003
        opt.sim_opt.p_21 = 4*opt.sim_opt.kappa_21/v_2; %kappa here is the membrane permeability
        v_1 = sqrt(4*opt.sim_opt.D_1/opt.sim_opt.sim_dt); %particle velocity in 1. This is adapted from Szafer A, Zhong J, Gore JC. 1995 and also Meier C, et al. 2003
        opt.sim_opt.p_12 = 4*opt.sim_opt.kappa_12/v_1; %kappa here is the membrane permeability
    end
    
    if strcmp(env, "one_gaussian_pool")
        opt.sim_opt.D_0 = 2e-9;
    end
    
    if strcmp(env, "three_gaussian_pools")
        opt.sim_opt.D_0 = [0.5e-9, 1e-9, 2e-9];
    end
    
    if strcmp(env, "one_circle")
        opt.sim_opt.one_cirlce_r = 5e-6;
        opt.sim_opt.D_0 = 2e-9;
    end
    
    
    if strcmp(env, "two_exchanging_pools")
        opt.sim_opt.f_1 = 0.7;
        opt.sim_opt.f_2 = 1-opt.sim_opt.f_1;
        opt.sim_opt.k = 20;
%         p_11 = @(t) opt.sim_opt.f_1 + opt.sim_opt.f_2*exp(-opt.sim_opt.k*t);
        p_12 = @(t) opt.sim_opt.f_2*(1 - exp(-opt.sim_opt.k*t));
%         p_22 = @(t) opt.sim_opt.f_2 + sim_opt.f_1*exp(-opt.sim_opt.k*t);
        p_21 = @(t) opt.sim_opt.f_1*(1 - exp(-opt.sim_opt.k*t));
        
        opt.sim_opt.p_12 = p_12(opt.sim_opt.sim_dt);
        opt.sim_opt.p_21 = p_21(opt.sim_opt.sim_dt);
        opt.sim_opt.D_1 = 1e-9;
        opt.sim_opt.D_2 = 2e-9;
    end

end