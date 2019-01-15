% Likelihood (up to normalizing constant) 
% of omega, rho, lambda(ii) given D

function log_likelihood = logLorl(omega,rho,lambda,obs_x,theta,...
    obs_y,sim_x,sim_t,sim_y)

%% For debugging: show inputs
% omega
% rho
% lambda
% obs_x
% theta
% obs_y

for ii = 1 : length(lambda)

    % Get m and n
    m = size(sim_x,1);
    n = size(obs_x,1);

    % Get observation covariance 
    if n>0
        error(['Current version of the code does not support '...
            'inclusion of field observations; instead, obs_x=[] required'])
        % Observation covariance should be specified
    else
        Sigma_y = [];
    end

    % Get simulation GP cov
    D_in_x = [ obs_x ; sim_x ];
    D_in_t = [ repmat(theta, n, 1) ; sim_t ];
    Sigma_eta = gp_cov(omega, D_in_x, D_in_x, rho, D_in_t, D_in_t,...
        lambda(ii),false);

    %% Get Sigma_D
    Sigma_D = Sigma_eta + padarray(Sigma_y,[m m],0,'post');
    % Add a nugget for computability
    WN = eye(size(Sigma_D)) * 10^(-4); Sigma_D = Sigma_D + WN;
    %rcond(Sigma_D)

    %% Get some needed values
    log_det_Sigma_D = logdet(Sigma_D);
    Sigma_D_inv = inv(Sigma_D);
    % DEPREC Sigma_D_inv_alt = Sigma_D\eye(size(Sigma_D));

    %% Get D
    D = [obs_y ; sim_y ] ;

    %% Get the likelihood
    omega_rho_lambda = [ omega rho lambda ] % Just to peek
    log_likelihood = -1/2 * log_det_Sigma_D - 1/2 * D' * Sigma_D_inv * D


end