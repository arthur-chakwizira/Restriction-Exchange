function resex_mc_one_circle_show_final_conf(R_new, r)
%show where particles were on last step

myfig;
x = linspace(-r,r);
c_p = real(sqrt( r^2 - (x).^2 )) ;
c_n =  real(-sqrt( r^2 - (x).^2 ));

fill(x, c_p, 'k', x, c_n ,'k')

plot(R_new(:, 1), R_new(:, 2), 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
end