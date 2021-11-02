function [R_new, states_new] = resex_mc_state_check_gaussian(R_old, states_old, sim_data)
%manage exchange process for two gaussian pools

    states_new = states_old;
    R_new = R_old;
    p_12 = sim_data.p_12;
    p_21 = sim_data.p_21;
    sigma_1 = sim_data.dr.sigma_1; %jump size intracellular
    sigma_2 = sim_data.dr.sigma_2; %jump size extracellular
    
    in_1 = (states_old == 1); %spins in pool 1
    move_or_not = rand(numel(in_1), 1); %transition probability
    can_move_to_2 = (in_1 & move_or_not < p_12); %move them to 2 if this is satisfied
    stay_in_1 = (in_1 & ~can_move_to_2); %otherwise keep them in 1
    
    R_new(can_move_to_2,:) = R_old(can_move_to_2, :) + randn(sum(can_move_to_2),2)*sigma_2; %update positions accordingly
    R_new(stay_in_1,:) = R_old(stay_in_1,:) + randn(sum(stay_in_1), 2)*sigma_1;
    
    %pool 2
    in_2 = (states_old == 2); %spins in pool 2
    move_or_not = rand(numel(in_2), 1); %%transition probability
    can_move_to_1 = (in_2 & move_or_not < p_21); %move them to 1 if this is satisfied
    stay_in_2 = (in_2 & ~can_move_to_1); %otherwise keep them in 2
    
    R_new(can_move_to_1, :) = R_old(can_move_to_1, :) + randn(sum(can_move_to_1), 2)*sigma_1;%update positions accordingly
    R_new(stay_in_2, :) = R_old(stay_in_2, :) + randn(sum(stay_in_2), 2)*sigma_2;
    
    %update states before next step
    states_new(stay_in_1) = 1; 
    states_new(can_move_to_2) = 2; 
    states_new(can_move_to_1) = 1; 
    states_new(stay_in_2) = 2;
    N_tot = sum(stay_in_1) + sum(can_move_to_2)+ sum(stay_in_2) + sum(can_move_to_1);
    if N_tot ~= numel(states_old); error('Conversation of particles violated.'); end
end