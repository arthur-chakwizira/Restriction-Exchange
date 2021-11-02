function [R_new, states_new] = resex_mc_state_check(sim_data, opt)
%This function accepts the previous spin positions and states and the new
%spin positions (preliminary) and re-calculates the new positions and states
%considering exchange between compartments.
%Input: exchange: struct with old positions, old states and new (preliminary)
%                 positions
%   :   restriction: struct with restriction-related parameters
%Output: R_new: final updated new positions (Nx2)
%      : states_new: final updated states (Nx1)

R_old_moved = sim_data.R_old_moved; %preliminary new positions
R_old = sim_data.R_old; %previous positions
states_old = sim_data.states_old; %previous states

p_12 = opt.sim_opt.p_12; %transmission probability intra- to extra
p_21 = opt.sim_opt.p_21; %transmission probability extra- to intra

r = opt.sim_opt.r; %cell radius
center_coords = sim_data.center_coords; %cell center positions

states_new = states_old; %initialise new states
R_new = zeros(size(R_old)); %initialise final positions

now_intra = false(size(R_old_moved,1), 1); %this array will be true for spins in intracellular compartment after update below

x_centres = center_coords.x;
y_centres = center_coords.y;
%go through all cells one at a time, with focus on intracellular
%compartment
for k = 1:size(x_centres, 1)
    for m = 1:size(y_centres, 1)
        c_x = x_centres(k, m);
        c_y = y_centres(k, m);
        in_cell = ( (R_old_moved(:,1) - c_x).^2 + (R_old_moved(:,2) - c_y).^2 ) < r^2; %point (x,y) is inside circle if (x-a)^2 + (y-b)^2 <= r^2
        now_intra(in_cell) = true; %these are particles now in this cell
        
        was_in_cell = ( (R_old(:,1) - c_x).^2 + (R_old(:,2) - c_y).^2 ) < r^2;
        was_intra = (states_old == 1);
        
        moved_legally = in_cell & was_in_cell; %these particles did not skip from one cell to another
        moved_illegally = in_cell & was_intra & ~was_in_cell; %these moved from one cell to another in time dt
        
        states_new(moved_legally|moved_illegally) = 1; %regardless of how they moved, they are intra now
        R_new(moved_legally, :) = R_old_moved(moved_legally, :); %accept current positions of legal immigrants
        R_new(moved_illegally, :) = R_old(moved_illegally, :); %restore illegal immigrants to previous positions
        
        
        %%
        was_extra = (in_cell & states_old == 2); %particles now in this cell and extracellular before must have come from the extracellular compartment
        allowed_into_1  = (rand(numel(was_extra),1) < p_21); %generate for each particle a random number in (0,1) and compare it to transition probability
        
        transmitted_into_1 = was_extra & allowed_into_1; %allow particles satisfying the condition above to enter this cell
        states_new(transmitted_into_1) = 1;%and change their states to intracellular
        
        R_new(transmitted_into_1, :) = R_old_moved(transmitted_into_1, :); %accept new positions, no need for complex rescaling
        
        rejected_from_1 = was_extra & ~allowed_into_1; %these particles are not allowed to complete transition into this cell
        states_new(rejected_from_1) = 2; %restore their states to  extracellular
        
        R_new(rejected_from_1, :) = R_old(rejected_from_1, :); %return particles to previous positions, skip compex rescaling
    end
end


now_extra = ~now_intra; %particles not found to be in any cell must be in the extracellular compartment
%go through all cells again, with focus on extracellular compartment
for k = 1:size(x_centres, 1)
    for m = 1:size(y_centres, 1)
        c_x = x_centres(k, m);
        c_y = y_centres(k, m);
        
        was_in_cell = ( (R_old(:,1) - c_x).^2 + (R_old(:,2) - c_y).^2 ) < r^2; %point (x,y) is inside circle if (x-a)^2 + (y-b)^2 <= r^2
        
        still_extra = (now_extra & states_old == 2); %these particles are extracellular both now and before
        states_new(still_extra) = 2; %keep their states extracellular
        R_new(still_extra, :) = R_old_moved(still_extra, :); %and accept their current positions
        
        was_intra = (was_in_cell & now_extra); %these particles were intracellular before and are now extracellular
        allowed_into_2  = (rand(numel(was_intra),1) < p_12);%generate for each particle a random number in (0,1) and compare it to transition probability
        
        transmitted_into_2 = was_intra & allowed_into_2; %allow particles satisfying the condition above to enter the extracellular compartment
        states_new(transmitted_into_2) = 2;%and change their states to extracellular
        R_new(transmitted_into_2, :) =  R_old_moved(transmitted_into_2, :);%simply accept the new positions
        
        rejected_from_2 = was_intra & ~allowed_into_2; %these particles not allowed to complete transition into extracellular
        states_new(rejected_from_2) = 1; %restore their states to intracellular;
        R_new(rejected_from_2, :) = R_old(rejected_from_2, :); %simply restore particle to previous positions
    end
end

end