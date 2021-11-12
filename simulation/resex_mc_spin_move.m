function [R_new, states_new, sim_data] = resex_mc_spin_move(sim_data, opt)
%This function accepts current positions and states and returns new
%positions and states after moving particles with fixed step lengths but 
%random directions. The new positions and states are returned after
%checking with the function "resex_mc_state_check".
%Input: sim_data: struct with old positions and states
%     : opt:  simulation options
%Output:  R_new: new positions (Nparticles x Ndim)
%      :  states_new: array of new states (Nparticles x 1)

R_old = sim_data.R_old; %previous positions
states_old = sim_data.states_old; %previous states

max_x = sim_data.max_dim.x; %maximum dimensions of substrate (boundaries)
max_y = sim_data.max_dim.y;
x_length = max_x(2) - max_x(1); %unit cell length in x
y_length = max_y(2) - max_y(1); %unit cell length in y


R_old_moved = R_old; %initialise new positions


%--RANDOM STEP LENGTHS, RANDOM DIRECTIONS------
R1 = R_old(states_old == 1, :); %previous intracellular positions
N1 = size(R1, 1);
R1_moved = R1 + randn(N1,2)*sim_data.sigma_1;

R2 = R_old(states_old == 2, :); %previous extracellular positions
N2 = size(R2, 1);
R2_moved = R2 + randn(N2,2)*sim_data.sigma_2;
%----------------------------------------------


%--FIXED STEP LENGTHS, RANDOM DIRECTIONS------
%{
% dr_1 = sim_data.dr.dr_1; %jump size intracellular (for fixed step in random directions)
% dr_2 = sim_data.dr.dr_2; %jump size extracellular (for fixed step in random directions)

R1 = R_old(states_old == 1, :); %previous intracellular positions
N1 = size(R1, 1);
x1_1 = 2*dr_1*rand(N1,1) - dr_1 +  R1(1:N1,1); %select random x-coordinate in [x0-dr_1, x0+dr_1]
rand_signs = randi(2, [N1,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
y1_1 = rand_signs.*sqrt(dr_1^2 - (x1_1 - R1(1:N1,1)).^2) + R1(1:N1,2);%calculate y-coordinate such that [x1_1,y1_1] lies on circle of radius dr_1
R1_moved = [x1_1, y1_1]; %then we have moved our spins in a random direction with a fixed jump length


R2 = R_old(states_old == 2, :); %previous extracellular positions
N2 = size(R2, 1);
x1_2 = 2*dr_2*rand(N2,1) - dr_2 +  R2(1:N2,1); %select random x-coordinate in [x0-dr_2, x0+dr_2]
rand_signs = randi(2, [N2,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
y1_2 = rand_signs.*sqrt(dr_2^2 - (x1_2 - R2(1:N2,1)).^2) + R2(1:N2,2);%calculate y-coordinate such that [x1_2,y1_2] lies on circle of radius dr_2
R2_moved = [x1_2, y1_2];%then we have moved our spins in a random direction with a fixed jump length
%}
%----------------------------------------------


R_old_moved(states_old == 1, :) = R1_moved; %update position array
R_old_moved(states_old == 2, :) = R2_moved;


sim_data.R_old_moved = R_old_moved; %pack new preliminary positions and send to state-checker function
[R_new, states_new] = resex_mc_state_check(sim_data, opt);


%now deal with particles that escaped the cell by imposing periodic boundaries

escaped_x_neg = R_new(:,1) < max_x(1); %left in negative x
R_new(escaped_x_neg, 1) = R_new(escaped_x_neg, 1) + x_length;
sim_data.h(escaped_x_neg, 1) = sim_data.h(escaped_x_neg, 1) - 1;

escaped_x_pos = R_new(:,1) >= max_x(2); %left in positive x
R_new(escaped_x_pos, 1) = R_new(escaped_x_pos, 1) - x_length;
sim_data.h(escaped_x_pos, 1) = sim_data.h(escaped_x_pos, 1) + 1;

escaped_y_neg = R_new(:,2) < max_y(1); %left in negative y
R_new(escaped_y_neg, 2) = R_new(escaped_y_neg, 2) + y_length;
sim_data.h(escaped_y_neg, 2) = sim_data.h(escaped_y_neg, 2) - 1;

escaped_y_pos = R_new(:,2) >= max_y(2); %left in positive y
R_new(escaped_y_pos, 2) = R_new(escaped_y_pos, 2) - y_length;
sim_data.h(escaped_y_pos, 2) = sim_data.h(escaped_y_pos, 2) + 1;

%ready for next time step
end