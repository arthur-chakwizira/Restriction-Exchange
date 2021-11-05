function     R_new = resex_mc_one_circle_restrict_particles(R_new, R_old, r)
%restrict particles to a circle of radius r

escaped =  ( (R_new(:,1)).^2 + (R_new(:,2)).^2  ) >= r^2;


R_new(escaped, :) = R_old(escaped, :);

end