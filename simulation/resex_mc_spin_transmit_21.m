function R_updated = resex_mc_spin_transmit_21(sim_data, opt)
%This function determines positions of particles transmitted from the
%extracellular compartment into the intracellular compartment
%Input: exchange struct with previous positions, current positions and
%       identities of admitted particles
%       restriction struct with cell position
%Output: R_updated (Nx2) containing transmitted positions

c_x = sim_data.c_x;
c_y = sim_data.c_y;

r = opt.sim_opt.r;

R_old = sim_data.R_old;
R_old_moved = sim_data.R_old_moved;
transmitted_into_1 = sim_data.transmitted_into_1;

D_1 = opt.sim_opt.D_1;
D_2 = opt.sim_opt.D_2;


if ~any(transmitted_into_1); R_updated = []; return; end

V0 = R_old(transmitted_into_1,:); %starting points
V1 = R_old_moved(transmitted_into_1,:); %end points
V0_1 = V1-V0; %trajectory vector from start to end

Vn = V0_1./(sqrt( V0_1(:,1).^2 + V0_1(:,2).^2)); %unit vectors along trajectory

h = linspace(0,1); %some parameter
Vx = V0(:,1) + h.*V0_1(:,1); %tracing particle x-positions fromm start to end
Vy = V0(:,2) + h.*V0_1(:,2); %same with y

d_center =  (Vx-c_x).^2 + (Vy-c_y).^2; %distances from the center squared
d_center_check = d_center < r^2;
Ventry = zeros(sum(transmitted_into_1),2); %will contain points of entry into intracellular
for row = 1:sum(transmitted_into_1)
    h_entry = h(find(d_center_check(row,:), 1, 'first')); %parameter value at entry point
    Ventry(row,:) = V0(row,:) + h_entry.*V0_1(row,:); %get entry point
end

V_in_1 = V1-Ventry; %vectors inside cell
d_in_1 = sqrt(V_in_1(:,1).^2 + V_in_1(:,2).^2); %path lengths in cell
d_scaled = d_in_1*sqrt(D_1/D_2); %scale path lengths with new diffusivity
V_final = Ventry+d_scaled.*Vn; %determine correct final positions in cell

R_updated = V_final; %update position vector
