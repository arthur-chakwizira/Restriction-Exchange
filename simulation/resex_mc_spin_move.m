function [R_new, states_new, sim_data] = resex_mc_spin_move(sim_data, opt)
%This function accepts current positions and states and returns new
%positions and states after moving particles with fixed step lengths but 
%random directions. The new positions and states are returned after
%checking with the function "res_exc_state_checker".
%Input: exchange: struct with old positions and states
%     : restriction:  struct with jump sizes
%Output:  R_new: new positions (Nx2)
%      :  states_new: array of new states (Nx1)

R_old = sim_data.R_old; %previous positions
states_old = sim_data.states_old; %previous states
sigma_1 = sim_data.dr.sigma_1; %jump size intracellular
sigma_2 = sim_data.dr.sigma_2; %jump size extracellular
max_x = sim_data.max_dim.x; %maximum dimensions of substrate (boundaries)
max_y = sim_data.max_dim.y;

x_length = max_x(2) - max_x(1);
y_length = max_y(2) - max_y(1);


R_old_moved = R_old; %initialise new positions

R1 = R_old(states_old == 1, :); %previous intracellular positions
N1 = size(R1, 1);
R1_moved = R1 + randn(N1,2)*sigma_1;

R2 = R_old(states_old == 2, :); %previous extracellular positions
N2 = size(R2, 1);
R2_moved = R2 + randn(N2,2)*sigma_2;

R_old_moved(states_old == 1, :) = R1_moved; %update position array
R_old_moved(states_old == 2, :) = R2_moved;


escaped_x_neg = R_old_moved(:,1) < max_x(1); %left in negative x
R_old_moved(escaped_x_neg, 1) = R_old_moved(escaped_x_neg, 1) + x_length;
escaped_x_pos = R_old_moved(:,1) >= max_x(2); %left in positive x
R_old_moved(escaped_x_pos, 1) = R_old_moved(escaped_x_pos, 1) - x_length;

escaped_y_neg = R_old_moved(:,2) < max_y(1); %left in negative y
R_old_moved(escaped_y_neg, 2) = R_old_moved(escaped_y_neg, 2) + y_length;
escaped_y_pos = R_old_moved(:,2) >= max_y(2); %left in positive y
R_old_moved(escaped_y_pos, 2) = R_old_moved(escaped_y_pos, 2) - y_length;


sim_data.R_old_moved = R_old_moved; %pack new preliminary positions and send to state-checker function

[R_new, states_new] = resex_mc_state_check(sim_data, opt);
end