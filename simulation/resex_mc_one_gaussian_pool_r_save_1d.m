function resex_mc_one_gaussian_pool_r_save_1d(r_fn, opt)
%This function simulates restricted diffusion with exchange on a 2D
%substrate. It can save mean spin positions every 1 millisecond.
%Key difference here is that time is now the third dimension. This should
%improve performance a lot

prog = 0;
fprintf(1,'Running simulation (one Gaussian pool): %3d%%\n',prog);

%% Retrieve parameters
T = opt.sim_opt.T; %total time
sim_dt = opt.sim_opt.sim_dt; %simulation time-step (less than 1 millisecond)
avg_dt = opt.sim_opt.avg_dt; %when to save mean positions;
N = round(opt.sim_opt.N_part); %intracellular population
D_0 = opt.sim_opt.D_0; %intracellular diffusivity
n_dim = 1; %number of dimensions, 1

sig_ma = sqrt(2*n_dim*D_0*sim_dt); %standard dev of steps in intracellular

R_old = randn(N, n_dim)*sqrt(2*D_0*sim_dt); % initial positions


%% Run simulation engine
sim_N_steps = round(T/sim_dt)+1; %number of steps in simulation

if opt.sim_opt.do_avg||opt.sim_opt.do_samp
    N_blocks = round(T/avg_dt)+1;%number of times average will be taken
    N_avg = round(avg_dt/sim_dt);%number of positions to average
    tmp_R = zeros(N, n_dim, N_avg); %storage of positions prior to averaging
    mean_R = zeros(N, n_dim, N_blocks); %storage of all average positions
    mean_R(:,:,1) = R_old; %initial values are positions at time 0
    counter = 1; %determines when it's time to average
    mean_counter = 1;%tracks number of averages taken
else
    mean_R =  zeros(N, n_dim, sim_N_steps);
    mean_R(:,:,1) = R_old; %initial values are positions at time 0
end

tic
loiter = opt.sim_opt.loiter; %allow some steps before starting to record
for step = 1:sim_N_steps+loiter
%     clc
%     disp(strcat('Step :', num2str(step), ' of :', num2str(sim_N_steps+loiter)))
    R_new = R_old + randn(N,n_dim)*sig_ma;
   
    R_old = R_new;
    
    
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
    prog = ( 100*(step/(sim_N_steps+loiter)) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end


    %% Save results and plot populations
    positions.opt = opt;
    positions.R = mean_R;
%     tic
    save(r_fn, 'positions','opt', '-v7.3')
     disp(" Alright, done.")
    toc
   
end
