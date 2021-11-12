function resex_mc_comptmt_vis(center_coords,r, max_dim)
%This function displays the 2D substrate (intra- and extracellular
%compartments), when evoked.
%Input: center_coords:  coordinates of centres of cells
%       r:  radius of cell
%       max_dim:  maximum extent of substrate
%Output: none

x_c = center_coords.x;
y_c = center_coords.y;

% 
% figure()
% axis equal
myfig;
grid off;
set(gcf, 'color', 'w');
set(gca,'box','on','linewidth',1,'layer','top', 'FontWeight', 'normal', 'FontSize', 12)
xlabel('distance [\mu m]')
ylabel('distance [\mu m]')
hold on
axis equal
for k = 1:size(x_c,1)
    for m = 1:size(x_c, 2)
        c_x = x_c(k, m)*1e6;
        c_y = y_c(k, m)*1e6;
        c = circle(c_x, c_y,r*1e6);
        fill(c(1,:), c(2,:), 'k', c(1,:), c(3,:) ,'k')
    end
end

% xmax = max_dim.x;
% ymax = max_dim.y;
% plot([xmax(1) xmax(1)], ymax, 'm-', xmax, [ymax(1) ymax(1)], 'm-', [xmax(2) xmax(2)], ymax, 'm-', xmax, [ymax(2) ymax(2)], 'm-')
% xlim(max_dim.x*1e6)
% ylim(max_dim.y*1e6)

    function c = circle(center_x, center_y, r) %circle generator
        x = linspace(center_x-r,center_x+r);
        c_p = real(sqrt( r^2 - (x-center_x).^2 )) + center_y;
        c_n =  real(-sqrt( r^2 - (x-center_x).^2 )) + center_y;
        c = [x; c_p; c_n];
    end

end