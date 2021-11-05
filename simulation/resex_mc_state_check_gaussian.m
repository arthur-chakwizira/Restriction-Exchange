function [R_new, states_new] = resex_mc_state_check_gaussian(R_old, states_old, sim_data)
%manage exchange process for two gaussian pools

    states_new = states_old;
    R_new = R_old;
    p_12 = sim_data.p_12;
    p_21 = sim_data.p_21;
    dr_1 = sim_data.dr.dr_1; %jump size intracellular
    dr_2 = sim_data.dr.dr_2; %jump size extracellular
    
    %pool 1
    in_1 = (states_old == 1); %spins in pool 1
    move_or_not = rand(numel(in_1), 1); %transition probability
    can_move_to_2 = (in_1 & move_or_not < p_12); %move them to 2 if this is satisfied
    stay_in_1 = (in_1 & ~can_move_to_2); %otherwise keep them in 1
    

N1 = sum(stay_in_1);
x1_1 = 2*dr_1*rand(N1,1) - dr_1 +  R_old(stay_in_1,1); %select random x-coordinate in [x0-dr_1, x0+dr_1]
rand_signs = randi(2, [N1,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
y1_1 = rand_signs.*sqrt(dr_1^2 - (x1_1 - R_old(stay_in_1,1)).^2) + R_old(stay_in_1,2);%calculate y-coordinate such that [x1_1,y1_1] lies on circle of radius dr_1
R_new(stay_in_1,:) = [x1_1, y1_1]; %then we have moved our spins in a random direction with a fixed jump length

N2 = sum(can_move_to_2);
x1_2 = 2*dr_2*rand(N2,1) - dr_2 +  R_old(can_move_to_2,1); %select random x-coordinate in [x0-dr_2, x0+dr_2]
rand_signs = randi(2, [N2,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
y1_2 = rand_signs.*sqrt(dr_2^2 - (x1_2 - R_old(can_move_to_2,1)).^2) + R_old(can_move_to_2,2);%calculate y-coordinate such that [x1_2,y1_2] lies on circle of radius dr_2
R_new(can_move_to_2, :) = [x1_2, y1_2];%then we have moved our spins in a random direction with a fixed jump length
%---------------------------------
    
    %pool 2
    in_2 = (states_old == 2); %spins in pool 2
    move_or_not = rand(numel(in_2), 1); %%transition probability
    can_move_to_1 = (in_2 & move_or_not < p_21); %move them to 1 if this is satisfied
    stay_in_2 = (in_2 & ~can_move_to_1); %otherwise keep them in 2
    
    
    N1 = sum(can_move_to_1);
x1_1 = 2*dr_1*rand(N1,1) - dr_1 +  R_old(can_move_to_1,1); %select random x-coordinate in [x0-dr_1, x0+dr_1]
rand_signs = randi(2, [N1,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
y1_1 = rand_signs.*sqrt(dr_1^2 - (x1_1 - R_old(can_move_to_1,1)).^2) + R_old(can_move_to_1,2);%calculate y-coordinate such that [x1_1,y1_1] lies on circle of radius dr_1
R_new(can_move_to_1,:) = [x1_1, y1_1]; %then we have moved our spins in a random direction with a fixed jump length

N2 = sum(stay_in_2);
x1_2 = 2*dr_2*rand(N2,1) - dr_2 +  R_old(stay_in_2,1); %select random x-coordinate in [x0-dr_2, x0+dr_2]
rand_signs = randi(2, [N2,1])-1; rand_signs(~rand_signs) = -1;%generate random signs to allow calculation of +/-sqrt(.)
y1_2 = rand_signs.*sqrt(dr_2^2 - (x1_2 - R_old(stay_in_2,1)).^2) + R_old(stay_in_2,2);%calculate y-coordinate such that [x1_2,y1_2] lies on circle of radius dr_2
R_new(stay_in_2, :) = [x1_2, y1_2];%then we have moved our spins in a random direction with a fixed jump length

    
    %update states before next step
    states_new(stay_in_1) = 1; 
    states_new(can_move_to_2) = 2; 
    states_new(can_move_to_1) = 1; 
    states_new(stay_in_2) = 2;
    N_tot = sum(stay_in_1) + sum(can_move_to_2)+ sum(stay_in_2) + sum(can_move_to_1);
    if N_tot ~= numel(states_old); error('Conversation of particles violated.'); end
end