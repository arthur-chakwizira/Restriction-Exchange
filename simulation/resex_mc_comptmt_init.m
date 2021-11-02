function sim_data = resex_mc_comptmt_init(opt)
%This function generates the initial distribution of particles and states
%in a 2D substrate modelling restricted diffusion and exchange.
%Input: exchange structure with exchange-related parameters
%     : restriction structure with restriction-related parameters
%Output:    R_start: initial positions of spins
%      :    center_coords: coordinates of centers of cells
%      :    max_dim: maximum extent of 2D substrate

r = opt.sim_opt.r;
show = opt.sim_opt.show;
N1 = round(opt.sim_opt.N_part*opt.sim_opt.f_1);
N2 = round(opt.sim_opt.N_part*opt.sim_opt.f_2);
N_cells = opt.sim_opt.N_cells;

center_x = -floor(N_cells/2)*2.5*r:2.5*r:floor(N_cells/2)*2.5*r; %gives cell separation of 0.5r
while mod((ceil(numel(center_x)/2)), 2) == 0
    N_cells = N_cells +1;
    center_x = -floor(N_cells/2)*2.5*r:2.5*r:floor(N_cells/2)*2.5*r; %gives cell separation of 0.5r
end
center_y = center_x;


center_coords = [center_x; center_y]; %positions of centers

cell_sep = 0.5*r;
center_shift = 1.5*r;

max_dim.x = [min(center_x)- r, max(center_x)+ r]; %tissue boundary in x-direction
max_dim.y = [min(center_y)- r, max(center_y) + center_shift + r]; %tissue boundary in y-direction

%intracellular
angles = rand(N1,1)*2*pi; %generate random angles in [0, 2*pi]
rs = rand(N1,1)*0.8*r; %and random radii in [0 0.8*r]
R_1 = [rs.*cos(angles), rs.*sin(angles)];  %distribute N1 spins uniformly on disk of radius 0.8r

%extracellular
angles = rand(N2,1)*2*pi;  %generate random angles in [0, 2*pi]
rs = rand(N2,1)*0.2*r + 1.2*r; %and random radii in [1.2r 1.4*r]
%  rs = rand(N2,1)*0.1*r + 1.1*r; %and random radii in [1.1r 1.2*r] %for spacing of 0.25r
 
R_2 = [rs.*cos(angles), rs.*sin(angles)];  %distribute N1 spins uniformly on disk of radius 0.8r

R_start = [R_1; R_2]; %initial positions


sim_data.R_old = R_start;
sim_data.center_coords = center_coords;
sim_data.max_dim = max_dim;
sim_data.cell_sep = cell_sep;
sim_data.center_shift = center_shift;


%visualise configuration if asked
if show
    resex_mc_comptmt_vis(center_coords, r, max_dim, sim_data);
    hold on
    set(gca, 'color', 'w');
    set(gca,'box','on','linewidth',1,'layer','top', 'FontWeight', 'normal', 'FontSize', 14)
    
    plot(R_start(1:N1,1), R_start(1:N1,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 1)
    plot(R_start(N1+1:N1+N2,1), R_start(N1+1:N1+N2,2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 1)
    
    xlabel('distance [m]')
    ylabel('distance [m]')
end



end