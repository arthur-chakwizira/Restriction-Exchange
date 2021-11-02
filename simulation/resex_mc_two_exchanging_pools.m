function resex_mc_two_exchanging_pools(r_fn, opt)
%This function simulates free diffusion with exchange in two infinite pools.
%It saves mean spin positions every 1 millisecond.
%Input: exchange parameter structure
%Output: None

%% Retrieve exchange parameters
T = exchange.T; %total time
sim_dt = exchange.sim_dt; %simulation time-step (less than 1 millisecond)
N1 = exchange.N1; %pool 1 population
N2 = exchange.N2; %pool 2 population
D_1 = exchange.D_1; %pool 1 diffusivity
D_2 = opt.D_2; %pool 2 diffusivity
N_steps = opt.N_steps;
n_dim = opt.n_dim; %number of dimensions, 2

%% Initialise pools
R = rand(N1+N2, n_dim)*sqrt(2*n_dim*D_1*sim_dt); %initial positions
dr_1 = sqrt(2*n_dim*D_1*sim_dt); %jump length in pool 1
dr_2 = sqrt(2*n_dim*D_2*sim_dt); %jump length in pool 2
p_12 = opt.p_12(sim_dt); %transition probability, from 1 to 2
p_21 = opt.p_21(sim_dt); %transition probability, from 2  to 1
states = ones(1,N1+N2); %initial states of spins, at time 0
states(N1+1:N1+N2) = states(N1+1:N1+N2)*2;

%% Run simulation engine
sim_N_steps = round(T/sim_dt)+1; %number of steps in simulation
N_blocks = round(opt.T/1e-3)+1; %number of times average will be taken
N_avg = round(opt.dt/opt.sim_dt);%number of positions to average
tmp_R = zeros(N_avg, N1+N2, 2); %storage of positions prior to averaging
mean_R = zeros(N_blocks, N1+N2,2);  %storage of all average positions
mean_R(1,:,:) = R; %initial values are positions at time 0
counter = 1; %determines when it's time to average
mean_counter = 1; %tracks number of averages taken

popn_1 = zeros(sim_N_steps,1);  %pool 1 population as function of time
popn_1(1) = N1;
popn_2 = zeros(sim_N_steps,1); %pool 2 population as function of time
popn_2(1) = N2;


loiter = 100;
for step = 1:sim_N_steps+loiter
    %pool 1
    in_1 = (states == 1); %spins in pool 1
    move_or_not = rand(1,numel(in_1)); %transition probability
    can_move_to_2 = (in_1 & move_or_not < p_12); %move them to 2 if this is satisfied
    stay_in_1 = (in_1 & ~can_move_to_2); %otherwise keep them in 1
    R(can_move_to_2,:) = exc_spin_mover(R(can_move_to_2,:), dr_2); %update positions accordingly
    R(stay_in_1,:) = exc_spin_mover(R(stay_in_1,:), dr_1);
    
    %pool 2
    in_2 = (states == 2); %spins in pool 2
    move_or_not = rand(1,numel(in_2)); %%transition probability
    can_move_to_1 = (in_2 & move_or_not < p_21); %move them to 1 if this is satisfied
    stay_in_2 = (in_2 & ~can_move_to_1); %otherwise keep them in 2
    R(can_move_to_1, :) = exc_spin_mover(R(can_move_to_1, :), dr_1); %update positions accordingly
    R(stay_in_2, :) = exc_spin_mover(R(stay_in_2, :), dr_2);
    
    %update states before next step
    states(stay_in_1) = 1;
    states(can_move_to_2) = 2;
    states(can_move_to_1) = 1;
    states(stay_in_2) = 2;
    N_tot = sum(stay_in_1) + sum(can_move_to_2)+ sum(stay_in_2) + sum(can_move_to_1);
    if N_tot ~= N1+N2; error('Conservation of particles violated.'); end
    
    if step > loiter %$$$$$$
        if opt.sim_opt.do_avg||opt.sim_opt.do_samp
            if counter <= N_avg %if it's not yet time to average
                tmp_R(:,:,counter) = R_new; %save positions to temporary array
                counter = counter + 1; %and keep counting
            else %if it's time to take average
                counter = 1; %reset counter
                mean_counter = mean_counter + 1; %update mean counter
                if opt.sim_opt.do_avg; tmp_mean_R = mean(tmp_R,3); end %compute temporal average of N_avg positions
                if opt.sim_opt.do_samp; tmp_mean_R = tmp_R(:, :, end); end %sample diffusion process every avg_dt
                mean_R(:,:,mean_counter) = tmp_mean_R; %save calculated mean
                tmp_R(:,:,counter) = R_new; %begin new batch of positions
                counter = 2; %update counter
            end
        else
            mean_R(:, :, step-loiter) = R_new;
        end
    end
    
    popn_1(step-loiter) = sum(states ==1); %update population arrays
    popn_2(step-loiter) = sum(states ==2);
end



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
toc
disp("Done.")
end