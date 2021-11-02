function resex_mc_final_r_plot(sim_data, opt)
%plots final positions of particles on substrate
    center_coords = sim_data.center_coords;
    r = opt.sim_opt.r;
    max_dim = sim_data.max_dim;
    R = sim_data.R_old;
    states = sim_data.states_old;
    
    resex_mc_comptmt_vis(center_coords, r, max_dim);
    hold on
    set(gca, 'color', 'w');
    set(gca,'box','on','linewidth',1,'layer','top', 'FontWeight', 'normal', 'FontSize', 14)
    
    plot(R(states == 1,1)*1e6, R(states==1,2)*1e6, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3)
    plot(R(states==2,1)*1e6, R(states==2,2)*1e6, 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 3)
    
    xlabel('distance [\mum]')
    ylabel('distance [\mum]')
    
    title('Final configuration')
end