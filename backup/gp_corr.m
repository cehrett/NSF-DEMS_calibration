% Calculate the correlation between (x,theta) and (x_s,theta_s) for the main
% Gaussian process of phase 1.

function corr_val = gp_corr(omega, rho, x, x_s, theta, theta_s)

corr_val = prod ( omega.^(4*(x-x_s).^2)) * prod (rho.^ (4*(theta-theta_s).^2));

end
