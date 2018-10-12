function results = MCMC_state_aware_discrep_true_fn ( settings ) 
% This function performs state-aware calibration while including a
% discrepancy function. The input settings is a structure which is itself
% output from another function.

%% Unpack settings
num_outputs = settings.num_outputs ; % Number of outputs from the model
lambda_delta = settings.lambda_delta ; % Marginal precision of discrepancy
x = settings.obs_locations ; % Control input points at which y is observed
xx = settings.obs_locns_with_dummy_vars_indexing_outputs; % x plus dum vars
y = settings.desired_observations ; % Target "observations" 
sigma2 = settings.observation_variance; % Known variance of observations
mu_theta = settings.mu_theta ; % mean of GP for state-aware theta1 
dim_theta1 = settings.dim_theta1; % Dimension of state-aware parameter
dim_theta2 = settings.dim_theta2; % Dim of non-state-aware parameter
dim_nu_theta = settings.dim_nu_theta; % Dim of inputs
dim_output = settings.dim_output; % Dim of output
a_theta = settings.lambda_theta_hypers(1); % Hyperparam, shape for gamma pr
b_theta = settings.lambda_theta_hypers(2); % Hyperparam, rate for gamma pr
Sigma_xi = settings.Sigma_xi; % Covariance for proposal dist for xi
Sigma_nu_theta = settings.Sigma_nu_theta; % Cov for prop dist for nu_theta
c_theta = settings.nu_theta_prior_param; % exp(-exp(nu_theta)) ~ Beta(1,c)
c_delta = settings.nu_delta_prior_param; % exp(-exp(nu_delta)) ~ Beta(1,c)

%% Create (known) observation variance matrix
Sigma = sigma2 * eye(length(y));

%% Create receptacles for sample draws
theta1_draws       = nan(M,numel(x)) ; 
xi_draws           = nan(M,dim_theta2) ;
nu_theta_draws     = nan(M,dim_nu_theta) ; 
lambda_theta_draws = nan(M,1) ;
nu_delta_draws     = nan(M,size(xx,2) * dim_output) ; 

%% Get initial draws of xi, nu_theta, lambda_theta, nu_delta
theta2_init             = rand(dim_theta2,1);
xi_draws(1,:)           = log(-log(theta2_init));
rho_theta_init          = betarnd(1,c_theta,dim_nu_theta,1);
nu_draws(1,:)           = log(-log(rho_theta_init));
lambda_theta_draws(1,:) = gamrnd(a_theta, 1/b_theta) ;
rho_delta_init          = betarnd(1,c_delta,size(xx,2)) ; 
nu_delta_draws(1,:)     = log(-log(rho_delta_init));

%% Get initial covariance matrices
% Although gp_cov accepts lambda as an input, here we give it lambda=1, so
% that we can handle lambda elsewhere in the calculations. (Giving gp_cov
% lambda=1 is effectively the same as removing lambda from gp_cov.)
R_delta = gp_cov(exp(-exp(nu_delta_draws(1,:))),xx,xx,0,0,0,1,false);
R_nu    = gp_cov(exp(-exp(nu_theta_draws(1,:))),x,x,0,0,0,1,false);

%% Get initial draw of theta1
theta1_draws(1,:) = mvnrnd(mu_theta*ones(length(x),1),...
    R_nu/lambda_theta_draws(1,:)) ;

%% Define scalar multipliers for adaptive proposals
mult_theta1       = 1 ;
mult_xi           = 1 ;
mult_nu_theta     = 1 ;
mult_lambda_theta = 1 ;
mult_nu_delta     = 1 ;

%% Define counts of accepted draws
% Will be used to monitor the MCMC and to tune the proposal distributions
accepted_theta1       = 0 ; 
accepted_xi           = 0 ;
accepted_nu_theta     = 0 ;
accepted_lambda_theta = 0 ;
accepted_nu_delta     = 0 ;

%% Get initial log likelihoods
% These will be updated whenever a proposed draw is accepted in MCMC.
logL_eta = logmvnpdf( y,...
    eta(xx,theta1_draws(1,:),exp(-exp(xi_draws(1,:)))),...
    Sigma + R_delta / lambda_delta) ;
logL_SA = logmvnpdf( link_fn(theta1_draws(1,:)),...
    mu_theta,...
    R_nu / lambda_theta_draws(1));
logL_prior_xi = sum(xi_draws(1,:)) - sum(exp(xi_draws(1,:)));
logL_prior_nu_theta = sum( (c_theta-1) * log(1-exp(-exp(nu_theta))));

%% MCMC
for ii = 2:M
    
    %% Draw theta1
    % Get proposal covariance (up to scale) for theta1. See Brown &
    % Atamturktur 2016 for details; they themselves are here following Neal
    % 1998.
    [U,D] = eig(R_nu) ;
    z = mvnrnd(0*ones(length(x),1),eye(length(x)));
    theta1_s = mult_theta1 * U * sqrt(D) * z + theta1_draws(ii-1,:) ;
    
    % Now find whether to accept or reject theta1_s. Each of the
    % log-likelihood of theta1 and theta1_s is found in two parts which are
    % summed; this is so that those parts can be used in later parts of the
    % MCMC rather than being recalculated.
    logL_theta1 = logL_eta + logL_SA; % log-likelihood of theta1
    
    logL_eta_s = logmvnpdf( y,...
        eta(xx,theta1_s,exp(-exp(xi_draws(ii-1,:)))),...
        Sigma + R_delta / lambda_delta);
    logL_SA_s = logmvnpdf( link_fn(theta1_draws(ii-1,:)),...
        mu_theta,...
        R_nu / lambda_theta_draws(ii-1));
    logL_theta1_s = logL_eta_s + logL_SA_s; % log-likelihood of theta1_s
    
    log_alpha = logL_theta1_s - logL_theta1 ; % log acceptance probability
    if log(rand) < log_alpha % if proposal is accepted
        accepted_theta1 = accepted_theta1 + 1 ;
        theta1_draws(ii,:) = theta1_s;
        logL_eta = logL_eta_a; % So that it can be used to draw xi
        logL_SA = logL_SA_s; % So that it can be used to draw nu_theta
    end
        
    %% Draw xi
    % Sample from proposal density:
    xi_s = mvnrnd(xi_draws(ii-1,:),mult_xi * Sigma_xi) ; 
    
    % Now find whether to accept or reject xi_s
    logL_xi = logL_eta + logL_prior_xi; % log-likelihood of xi
    
    logL_eta_s = logmvnpdf( y,...
        eta(xx,theta1_draws(ii,:),exp(-exp(xi_s))),...
        Sigma + R_delta / lambda_delta) ;
    logL_prior_xi_s = sum(xi_s) - sum(exp(xi_s));
    logL_xi_s = logL_eta_s + logL_prior_xi_s; % log-likelihood of xi_s
    
    log_alpha = logL_xi_s - logL_xi ; % log acceptance probability
    if log(rand) < log_alpha % if proposal is accepted
        accepted_xi = accepted_xi + 1 ;
        xi_draws(ii,:) = xi_s ; 
        logL_eta = logL_eta_s ;
        logL_prior_xi = logL_prior_xi_s ;
    end
    
    %% Draw nu_theta
    % Sample from proposal density and make corresponding covariance matrix
    nu_theta_s = ...
        mvnrnd(nu_theta_draws(ii-1,:), mult_nu_theta * Sigma_nu_theta);
    R_nu_s = gp_cov(exp(-exp(nu_theta_s)),x,x,0,0,0,1,false);
    
    % Now find whether to accept or reject nu_theta_s
    logL_nu_theta = logL_SA + logL_prior_nu_theta ; % log-L of nu_theta
    
    logL_SA_s = logmvnpdf( link_fn(theta1_draws(ii,:)),...
        mu_theta,...
        R_nu_s / lambda_theta_draws(ii-1));
    logL_prior_nu_theta_s = sum( (c_theta-1) * log(1-exp(-exp(nu_theta))));
    logL_nu_theta_s = logL_SA_s +logL_prior_nu_theta_s;%log-L of nu_theta_s
    
    log_alpha = logL_nu_theta_s - logL_nu_theta ; % log acceptance prob.
    if log(rand) < log_alpha % if proposal is accepted
        accepted_nu_theta = accepted_nu_theta + 1 ;
        nu_theta_draws(ii,:) = nu_theta_s ; 
        logL_SA = logL_SA_s ; 
        logL_prior_nu_theta = logL_prior_nu_theta_s ; 
        R_nu = R_nu_s ; 
    end
    
    %% Draw lambda_theta
    % Sample from proposal density
    
    
end


end