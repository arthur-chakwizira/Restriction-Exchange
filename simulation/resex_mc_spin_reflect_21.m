function R_updated = resex_mc_spin_reflect_21(sim_data, opt)
%This function determines positions of particles reflected from the
%extracellular compartment.
%Input: exchange struct with previous positions, current positions and
%       identities of rejected particles
%       restriction struct with cell position
%Output: R_updated (Nx2) containing reflected positions

c_x = sim_data.c_x;
c_y = sim_data.c_y;

r = opt.sim_opt.r;

R_old = sim_data.R_old;
R_old_moved = sim_data.R_old_moved;
rejected_from_2 = sim_data.rejected_from_2;

if ~any(rejected_from_2); R_updated = []; return; end

V0 = R_old(rejected_from_2,:); %starting points
V1 = R_old_moved(rejected_from_2,:); %end points
V0_1 = V1-V0; %trajectory vector from start to end

Vn = V0_1./(sqrt( V0_1(:,1).^2 + V0_1(:,2).^2)); %unit vectors along trajectory
h = linspace(0,1); %some parameter
Vx = V0(:,1) + h.*V0_1(:,1); %tracing particle x-positions from start to end
Vy = V0(:,2) + h.*V0_1(:,2); %same with y
d_center =  (Vx-c_x).^2 + (Vy-c_y).^2; %distances from the center, squared
d_center_check = (d_center < r^2);
Ventry = zeros(sum(rejected_from_2), 2); %will contain points of entry into extracellular
for row = 1:sum(rejected_from_2)
    h_entry = h(find(d_center_check(row,:), 1, 'last')); %parameter value at entry point
    Ventry(row,:) = V0(row,:) + h_entry*V0_1(row,:); %get entry points
end
V_in_2 = V1-Ventry; %vectors inside extracellular
d_in_2 = sqrt(V_in_2(:,1).^2 + V_in_2(:,2).^2); %path lengths in extracellular
V_final = Ventry - d_in_2.*Vn; %determine correct final positions inside cell

failed = (V_final(:,1)-c_x).^2 + (V_final(:, 2)-c_y).^2 >= r^2; %happens if exit trajectory is tangential to boundary
V_final(failed,:) = R_old(failed,:);
R_updated = V_final; %update position vector
end


