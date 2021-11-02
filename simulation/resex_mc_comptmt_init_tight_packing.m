function sim_data = resex_mc_comptmt_init_tight_packing(opt)
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


N_cells = 9;
max_dim.x = [0, 0]; %unit cell boundary in x-direction
max_dim.y = [1, 1]; %unit cell boundary in y-direction
cell_sep = 2.5e-6;
if cell_sep < sqrt(2*opt.sim_opt.n_dim*opt.sim_opt.D)
center_shift = r + cell_sep;
cell_room = 2*r + cell_sep; %space taken up by one cell

while ~isequal(max_dim.x, max_dim.y)
center_x = -floor(N_cells/2)*cell_room:cell_room:floor(N_cells/2)*cell_room; %gives cell separation of cell_sep

N_cells = N_cells+1;
if mod((ceil(numel(center_x)/2)), 2) == 0; continue; end

center_y = center_x;


%make substate grid
center_coords.x = repmat(center_x, numel(center_x),1);
center_coords.y = repmat(center_y, numel(center_y),1)';


% shift every second column for tighter packing
for k = 1:numel(center_x)
        if mod(k,2) == 0;center_coords.y(:,k) = center_coords.y(:,k) + center_shift; end
end



%slice grid into unit cell 
origin_column = center_coords.y(:, ceil(N_cells/2));

garbage = abs(origin_column) > 6*r + 3*cell_sep; %remove these rows; we want a 6x6 cell box



center_coords.y(garbage,:) = [];
center_coords.x(garbage,:) = [];

garbage = (abs(center_coords.y) > 6*r + 3*cell_sep);
center_coords.y(garbage) = NaN;
center_coords.x(garbage) = NaN;

any_row = center_coords.x(1,:);
garbage = abs(any_row) > 6*r + 3*cell_sep;
center_coords.x(:, garbage) = [];
center_coords.y(:, garbage) = [];

max_dim.x = [min(center_coords.x(:)), max(center_coords.x(:))]; %unit cell boundary in x-direction
max_dim.y = [min(center_coords.y(:)), max(center_coords.y(:))]; %unit cell boundary in y-direction

end
%okay the slicing here is actually quite cumbersome. Might be worth
%considering a differnt packing geometry


%now distribute particles all over substrate, not just in the centre

%%
R_1 = zeros(N1, 2);
R_2 = zeros(N2,2);

c1 = 0;
c2 = 0;
while (c1 < N1)
    x = max_dim.x(1) + rand*(max_dim.x(2)-max_dim.x(1));
    y = max_dim.y(1) + rand*(max_dim.y(2)-max_dim.y(1));
    if any(( (x - center_coords.x).^2 + (y - center_coords.y).^2 ) < r^2, 'all') %coordinates are in some cell
        c1 = c1+1;
        R_1(c1, :) = [x, y];
    end
end


while (c2 < N2)
    x = max_dim.x(1) + rand*(max_dim.x(2)-max_dim.x(1));
    y = max_dim.y(1) + rand*(max_dim.y(2)-max_dim.y(1));
    if ~any(( (x - center_coords.x).^2 + (y - center_coords.y).^2 ) < r^2, 'all') %coordinates outside any cell
        c2 = c2+1;
        R_2(c2, :) = [x, y];
    end
end


R_start = [R_1; R_2]; %initial positions
%%

%visualise configuration if asked
sim_data.R_old = R_start;
sim_data.center_coords = center_coords;
sim_data.max_dim = max_dim;

if show
    resex_mc_comptmt_vis(center_coords, r, max_dim);
    hold on
    set(gca, 'color', 'w');
    set(gca,'box','on','linewidth',1,'layer','top', 'FontWeight', 'normal', 'FontSize', 14)
    
    plot(R_start(1:N1,1)*1e6, R_start(1:N1,2)*1e6, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 1)
    plot(R_start(N1+1:N1+N2,1)*1e6, R_start(N1+1:N1+N2,2)*1e6, 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 1)
end


end