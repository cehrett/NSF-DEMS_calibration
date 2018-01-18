%% Covariance parameter settings
omega = [0.2 0.82];
rho   = [0.39 0.99];
lambda = 0.5; 

x = linspace(0,1,3);
R = gp_cov(omega,x,x,rho,x,x,lambda);

save('test_R','R');
