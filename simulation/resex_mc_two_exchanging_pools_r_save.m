function resex_mc_two_exchanging_pools_r_save(r_fn, opt)
%This function simulates restricted diffusion with exchange on a 2D
%substrate. It can save mean spin positions every 1 millisecond.
tic
prog = 0;
fprintf(1,'Running simulation (two exchanging Gaussian pools): %3d%%\n',prog);
%% Retrieve parameters
T = opt.sim_opt.T; %total time
sim_dt = opt.sim_opt.sim_dt; %simulation time-step (less than 1 millisecond)
avg_dt = opt.sim_opt.avg_dt; %when to save mean positions;
N1 = round(opt.sim_opt.f_1*opt.sim_opt.N_part); %intracellular population
N2 = round(opt.sim_opt.f_2*opt.sim_opt.N_part); %extracellular population
D_1 = opt.sim_opt.D_1; %intracellular diffusivity
D_2 = opt.sim_opt.D_2; %extracellular diffusivity
% n_dim = opt.sim_opt.n_dim; %number of dimensions, 2

%% Initialise
states_old = ones(N1+N2,1);  %initial states of intracellular spins (at time 0)
states_old(N1+1:N1+N2) = states_old(N1+1:N1+N2)*2; %initial states of extracellular spins
R_old = randn(N1+N2, 2);
sigma_1 = sqrt(2*D_1*sim_dt); %standard dev of steps in intracellular
sigma_2 = sqrt(2*D_2*sim_dt); %standard dev of steps in extracellular
dr.sigma_1 = sigma_1; %pack parameters for simulation engine
dr.sigma_2 = sigma_2;

sim_data.dr = dr;
sim_data.p_12 = opt.sim_opt.p_12;
sim_data.p_21 = opt.sim_opt.p_21;


%% Run simulation engine
sim_N_steps = round(T/sim_dt)+1; %number of steps in simulation

if opt.sim_opt.do_avg||opt.sim_opt.do_samp
    N_blocks = round(T/avg_dt)+1;%number of times average will be taken
    N_avg = round(avg_dt/sim_dt);%number of positions to average
    tmp_R = zeros(N1+N2, 2, N_avg); %storage of positions prior to averaging
    mean_R = zeros(N1+N2,2, N_blocks); %storage of all average positions
    mean_R(:,:,1) = R_old; %initial values are positions at time 0
    counter = 1; %determines when it's time to average
    mean_counter = 1;%tracks number of averages taken
else
    mean_R =  zeros(N1+N2, 2, sim_N_steps);
    mean_R(:,:,1) = R_old; %initial values are positions at time 0
end

popn_1 = zeros(sim_N_steps,1); %intracellular population as function of time
popn_1(1) = N1;
popn_2 = zeros(sim_N_steps,1); %extracellular population as function of time
popn_2(1) = N2;

loiter = opt.sim_opt.loiter; %allow some steps before starting to record
for step = 1:sim_N_steps+loiter
        
    [R_new, states_new] = resex_mc_state_check_gaussian(R_old, states_old, sim_data); 
    R_old = R_new;
    states_old = states_new;
    
    
    if step > loiter %$$$$$$
        if opt.sim_opt.do_avg||opt.sim_opt.do_samp
            if counter <= N_avg %if it's not yet time to average
                tmp_R(:,:, counter) = R_new; %save positions to temporary array
                counter = counter + 1; %and keep counting
            else %if it's time to take average
                counter = 1; %reset counter
                mean_counter = mean_counter + 1; %update mean counter
                if opt.sim_opt.do_avg; tmp_mean_R = mean(tmp_R,3); end %compute average of N_avg positions
                if opt.sim_opt.do_samp; tmp_mean_R = tmp_R(:, :, end); end %sample diffusion process every avg_dt
                mean_R(:,:, mean_counter) = tmp_mean_R; %save calculated mean
                tmp_R(:,:,counter) = R_new; %begin new batch of positions
                counter = 2; %update counter
            end
        else
            mean_R(:, :, step-loiter) = R_new;
        end
        
        popn_1(step-loiter) = sum(states_new ==1);%update population arrays
        popn_2(step-loiter) = sum(states_new ==2);
        
    end
    
    
        
    prog = ( 100*(step/(sim_N_steps+loiter)) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end %$$$$$$

try close(writerObj); catch; end

%% Save results and plot populations
positions.opt = opt;
% positions.sim_data = sim_data;
positions.R = mean_R;
positions.popn_1 = popn_1;
positions.popn_2 = popn_2;

save(r_fn, 'positions', 'opt','-v7.3')


time = linspace(0,T,sim_N_steps);
% figure()
% plot(time, popn_1, 'r-', 'LineWidth', 2)
% hold on
% plot(time, popn_2, 'c-', 'LineWidth', 2)
% legend({'N1', 'N2'})
% title('Population visualisation')
%
% popn_1 = popn_2;
save('resex_mc_popn_1_gaussian.mat', 'popn_1', 'popn_2', 'time')

disp(" Done")
toc
end
