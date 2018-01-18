% Likelihood (up to normalizing constant) 
% of omega, rho, lambda(ii) given D


function likelihood = Lorl(omega,rho,lambda,obs_x,theta,...
    obs_y,sim_x,sim_t,sim_y,c)

%% For debugging: show inputs
% omega
% rho
% lambda
% obs_x
% theta
% obs_y

if c == 0 % c is a scalar that is used for numerical stability,
    % by multiplying the covariance matrix prior to finding its det
    c = 10^3.6;
end

for ii = 1 : length(lambda)

% Get m and n
m = size(sim_x,1);
n = size(obs_x,1);

% Get observation covariance 
if n>0
    Sigma_y = cov(obs_y);
else
    Sigma_y = [];
end

% Get simulation GP cov
D_in_x = [ obs_x ; sim_x ];
D_in_t = [ repmat(theta, n, 1) ; sim_t ];
Sigma_eta = gp_cov(omega, D_in_x, D_in_x, rho, D_in_t, D_in_t, lambda(ii));

%% Get Sigma_D
Sigma_D = Sigma_eta + padarray(Sigma_y,[m m],0,'post');
% Add a nugget for computability
WN = eye(size(Sigma_D)) * 10^(-4); Sigma_D = Sigma_D + WN;
%rcond(Sigma_D)

%% Get some needed values
det_Sigma_D = det(c * Sigma_D); % c is used for numerical stability
% DEPREC Use LU decomp to stay on log scale
% [L,U] = lu(Sigma_D);
% log_det_Sigma_D = sum(log(diag(U)));
Sigma_D_inv = inv(Sigma_D);
% DEPREC Sigma_D_inv_alt = Sigma_D\eye(size(Sigma_D));

%% Get D
D = [obs_y ; sim_y ] ;

%% Get the likelihood
omega_rho_lambda = [ omega rho lambda ]
likelihood = det_Sigma_D^(-1/2) * exp(-1/2 * D' * Sigma_D_inv * D)
% likelihood = mvnpdf(D,0,Sigma_D)
% DEPREC Use log scale
% loglik(ii) = -1/2 * log_det_Sigma_D - 1/2 * D' * Sigma_D_inv * D;
%loglik_alt = -1/2 * log_det_Sigma_D - 1/2 * D' * Sigma_D_inv_alt * D


end