function [center_coords, max_dim] = resex_mc_comptmt_grow(sim_data, opt)
r = opt.sim_opt.r;
N_cells = opt.sim_opt.N_cells;

accept = false; %dont just grow substrate willy-nilly
while ~accept
    if mod(N_cells, 2) ~= 0 && (mod(ceil(N_cells/2), 2)) ~= 0
        accept = true;
    else
        N_cells = N_cells +1;
    end
end

center_shift = sim_data.center_shift;

center_x = -floor(N_cells/2)*2.5*r:2.5*r:floor(N_cells/2)*2.5*r; %gives cell separation of 0.5r
center_y = center_x;
center_coords = [center_x; center_y]; %positions of centers

max_dim.x = [min(center_x)- r, max(center_x)+ r]; %tissue boundary in x-direction
max_dim.y = [min(center_y)- r, max(center_y) + center_shift + r]; %tissue boundary in y-direction


end