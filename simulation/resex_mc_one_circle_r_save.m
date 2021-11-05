function resex_mc_one_circle_r_save(r_fn, opt)
%This function simulates restricted diffusion with exchange on a circle
%. It can save mean spin positions every 1 millisecond.
tic
prog = 0;
fprintf(1,'Running simulation (one circle): %3d%%\n',prog);
%% Retrieve parameters
T = opt.sim_opt.T; %total time
sim_dt = opt.sim_opt.sim_dt; %simulation time-step (less than 1 millisecond)
avg_dt = opt.sim_opt.avg_dt; %when to save mean positions;
N = round(opt.sim_opt.N_part); %intracellular population
D_0 = opt.sim_opt.D_0; %intracellular diffusivity
n_dim = opt.sim_opt.n_dim; %number of dimensions, 2

dr = sqrt(2*n_dim*D_0*sim_dt); %standard dev of steps in intracellular

r = opt.sim_opt.one_cirlce_r; %radius of the circle

%place particles initially in cirlce
X_0 = -r + 2*r*rand(N, 1);
Y_0 = sqrt((r)^2- X_0.^2);
Y_0 = -Y_0 + 2*rand(size(Y_0)).*Y_0;
R_old = [X_0 Y_0]*0; % initial positions

% plot(X_0, Y_0, 'k+')
% error('Sir')
%% Run simulation engine
sim_N_steps = round(T/sim_dt)+1; %number of steps in simulation

if opt.sim_opt.do_avg||opt.sim_opt.do_samp
    N_blocks = round(T/avg_dt)+1;%number of times average will be taken
    N_avg = round(avg_dt/sim_dt);%number of positions to average
    tmp_R = zeros(N, 2, N_avg); %storage of positions prior to averaging
    mean_R = zeros(N,2, N_blocks); %storage of all average positions
    mean_R(:,:, 1) = R_old; %initial values are positions at time 0
    counter = 1; %determines when it's time to average
    mean_counter = 1;%tracks number of averages taken
else
    mean_R =  zeros(N, 2, sim_N_steps);
    mean_R(:,:, 1) = R_old; %initial values are positions at time 0
end

loiter = opt.sim_opt.loiter; %allow some steps before starting to record
for step = 1:sim_N_steps+loiter
%     clc
%     disp(strcat('Step :', num2str(step), ' of :', num2str(sim_N_steps+loiter)))
%     R_new = R_old + randn(N,2)*sig_ma;
    
    x1 = 2*dr*rand(N,1) - dr +  R_old(:,1); %select random x-coordinate in [x0-dr, x0+dr]
    rand_signs = randi(2, [N,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
    y1 = rand_signs.*sqrt(dr^2 - (x1 - R_old(:,1)).^2) + R_old(:,2);%calculate y-coordinate such that [x1,y1] lies on circle of radius dr
    R_new = [x1, y1]; %then we have moved our spins in a random direction with a fixed jump length

    
    %restrict particles to circle
    R_new = resex_mc_one_circle_restrict_particles(R_new, R_old, r);
    
    R_old = R_new;
    
        
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
                tmp_R(:,:, counter) = R_new; %begin new batch of positions
                counter = 2; %update counter
            end
        else
            mean_R(:, :, step-loiter) = R_new;
        end
    end
    prog = ( 100*(step/(sim_N_steps+loiter)) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end


    %% Save results and plot populations
    positions.opt = opt;
    positions.R = mean_R;
    
    resex_mc_one_circle_show_final_conf(R_new, r)
    
    save(r_fn, 'positions','opt', '-v7.3')
     disp(" Alright, done.")
    toc
end
