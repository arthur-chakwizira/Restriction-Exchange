function n_cells = resex_mc_get_num_cells(dt, T, D, d, n_dim, loiter)
%compute number of cells required for simulation
%conservative (it overestimates the number of cells to be safe)

t_tot = T + loiter*dt;
max_range = sqrt(2*n_dim*D*(t_tot));

n_cells = ceil(3*max_range/d);

while mod(ceil(n_cells/2), 2) == 0
    n_cells = n_cells + 1;
end
end