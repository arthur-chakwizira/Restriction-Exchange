function resex_mc_r_save(r_fn, opt)
%This function simulates restricted diffusion with exchange on a 2D
%substrate.
tic

%% Retrieve parameters
T = opt.sim_opt.T; %total time
sim_dt = opt.sim_opt.sim_dt; %simulation time-step (less than 1 millisecond)
avg_dt = opt.sim_opt.avg_dt; %when to save mean positions;
N1 = round(opt.sim_opt.f_1*opt.sim_opt.N_part); %intracellular population
N2 = round(opt.sim_opt.f_2*opt.sim_opt.N_part); %extracellular population
D_1 = opt.sim_opt.D_1; %intracellular diffusivity
D_2 = opt.sim_opt.D_2; %extracellular diffusivity
n_dim = opt.sim_opt.n_dim; %number of dimensions, 2
r = opt.sim_opt.r; %compartment radius

%% Initialise compartments
% sim_data = resex_mc_comptmt_init(opt); %generate initial positions and get substrate extent
sim_data = resex_mc_comptmt_init_tight_packing(opt);
states_old = ones(N1+N2,1);  %initial states of intracellular spins (at time 0)
states_old(N1+1:N1+N2) = states_old(N1+1:N1+N2)*2; %initial states of extracellular spins
sim_data.states_old = states_old;
sim_data.sigma_1 = sqrt(2*1*D_1*sim_dt);
sim_data.sigma_2 = sqrt(2*1*D_2*sim_dt);
dr_1 = sqrt(2*n_dim*D_1*sim_dt); %jump length in intracellular
dr_2 = sqrt(2*n_dim*D_2*sim_dt); %jump length in extracellular
dr.dr_1 = dr_1; %pack parameters for simulation engine
dr.dr_2 = dr_2;
sim_data.dr = dr;


%% Prepare for simulation visualisation
film = opt.sim_opt.film;
if film.make
    close all;  resex_mc_comptmt_vis(sim_data.center_coords,r, sim_data.max_dim); hold on;
    film.plot_1 = plot( NaN,  NaN,'ro','MarkerFaceColor','r', 'MarkerSize', 6);
    film.plot_2 = plot( NaN,  NaN,'co','MarkerFaceColor','c', 'MarkerSize', 6);
    title(["PARTICLE DYNAMICS"; " "], 'Color', 'r')
    if film.save_film
        writerObj = VideoWriter(film.film_title);
        writerObj.FrameRate = film.frame_rate;
        film.writerObj = writerObj;
        open(writerObj);
    end
end

%% Run simulation engine
sim_N_steps = round(T/sim_dt)+1; %number of steps in simulation

R = zeros(opt.sim_opt.N_part, opt.sim_opt.n_dim, sim_N_steps); %will contain entire position history
states = zeros(opt.sim_opt.N_part,sim_N_steps); %will contain entire compartment history

popn_1 = zeros(sim_N_steps,1); %intracellular population as function of time
popn_1(1) = N1;
popn_2 = zeros(sim_N_steps,1); %extracellular population as function of time
popn_2(1) = N2;

h = zeros(opt.sim_opt.N_part, 2); %keep track of jumps
sim_data.h = h;

loiter = opt.sim_opt.loiter; %allow some steps before starting to record

% allow delay before recording
fprintf(1, 'Preparing for simulation (two compartments)... \n')
for step = 1:loiter
    [R_new, states_new, sim_data] = resex_mc_spin_move(sim_data, opt); %move spins
    sim_data.R_old = R_new;
    sim_data.states_old = states_new;
end


% now record
fprintf(1,'Running simulation : %3d%%\n',0);
for step = 1:sim_N_steps
    [R_new, states_new, sim_data] = resex_mc_spin_move(sim_data, opt); %move spins
    sim_data.R_old = R_new;
    sim_data.states_old = states_new;
    
    R(:,:,step) = R_new + sim_data.h*sim_data.cell_length;
    states(:, step) = states_new;
    %______________________________________________________________________
    if film.make
        film.step = step; film.R_new = R_new; film.states_new = states_new; film.last_step = sim_N_steps;
        resex_mc_spin_visualiser(film)
    end
    %______________________________________________________________________
    popn_1(step) = sum(states_new ==1);%update population arrays
    popn_2(step) = sum(states_new ==2);

    fprintf(1,'\b\b\b\b%3.0f%%', 100*(step/(sim_N_steps)) );
end
try close(writerObj); catch; end


%% sample/average if needed
fprintf(1,'\nFinishing up... \n');

if opt.sim_opt.do_avg||opt.sim_opt.do_samp
    N_blocks = round(T/avg_dt);%number of times average will be taken
    N_avg = round(avg_dt/sim_dt);%number of positions to average
    tmp_R = zeros(N1+N2, 2, N_avg); %storage of positions prior to averaging
    mean_R = zeros(N1+N2,2, N_blocks); %storage of all average positions
%     mean_R(:,:,1) = sim_data.R_old; %initial values are positions at time 0
    counter = 1; %determines when it's time to average
    mean_counter = 0;%tracks number of averages taken
    
    for step = 1:sim_N_steps
        if counter <= N_avg %if it's not yet time to average
            tmp_R(:,:, counter) = R(:,:,step); %save positions to temporary array
            counter = counter + 1; %and keep counting
        else %if it's time to take average
            counter = 1; %reset counter
            mean_counter = mean_counter + 1; %update mean counter
            if opt.sim_opt.do_avg; tmp_mean_R = mean(tmp_R,3); end %compute average of N_avg positions
            if opt.sim_opt.do_samp; tmp_mean_R = tmp_R(:, :, end); end %sample diffusion process every avg_dt
            mean_R(:,:, mean_counter) = tmp_mean_R; %save calculated mean
            tmp_R(:,:,counter) = R(:,:,step); %begin new batch of positions
            counter = 2; %update counter
        end
    end
else
    mean_R = R;
end


%% Save results and plot populations
positions.opt = opt;
positions.R = mean_R;
positions.popn_1 = popn_1;
positions.popn_2 = popn_2;
positions.sim_data = sim_data;

save(r_fn, 'positions','-v7.3')


time = linspace(0,T,sim_N_steps);
figure('Color', 'w')
plot(time, popn_1, 'r-', 'LineWidth', 2)
hold on
plot(time, popn_2, 'c-', 'LineWidth', 2)
legend({'N1', 'N2'})
title('POPULATION VISUALISATION')
set(gca, 'fontsize', 20)
grid on
xlabel('time [s]')

save('resex_mc_popn_1.mat', 'popn_1', 'popn_2', 'time')

%show where particles ended up on grid
sim_data.R_old = R(:, :, end);
resex_mc_final_r_plot(sim_data, opt)
fprintf(1,'Done \n');
toc
end

