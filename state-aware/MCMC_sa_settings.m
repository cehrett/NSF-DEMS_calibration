function settings = MCMC_sa_settings (...
    desired_obs,lambda_delta,cntrl_input,varargin)
% This function packages settings for the state-aware MCMC calibration
% using MCMC_state_aware_discrep_true_fn().
% Required inputs: desired_obs is the desired observation, lambda_delta is
% the estimate of the distance of the desired observation from the Pareto
% front, cntrl_input is the vector or matrix of control input settings, 
% on the normalized scale, at which the desired observation is to be 
% "observed" and at which values of theta1 are to be estimated. Further 
% parameters may be defined as follows.
% 
% Parameter/Value pairs:
% 'M':
% Total number of iterations, including burn-in. Default: 2e4.
%
% 'burn_in':
% Total number of iterations to treat as burn-in. Default: 1/2 * M.
% 
% 'nugget':
% Size of the irreducible error variance. Default: 0.05 (standardized
% scale).
%
% 'mu_theta': 
% Mean of GP for state-aware calibration. Default: 0.5.
% 
% 'dim_theta1':
% Dimension of the state-aware calibration parameter. Default: 1.
% 
% 'dim_theta2':
% Dimension of the non-state-aware calibration parameter. Default: 1.
% 
% 'lambda_theta_hypers':
% Length-2 vector of gamma hyperparameters for the prior distribution on
% lambda_theta, using shape-rate parameterization. Default: [0.1, 0.1].
%
% 'Sigma_xi':
% Covariance for initial proposal dist for xi. Early in the burn-in period,
% this will be replaced with an adaptive covariance matrix. Default:
% identity matrix.
%
% 'Sigma_nu_theta':
% Covariance for proposal distribution for nu_theta. Early in the burn-in 
% period, this will be replaced with an adaptive covariance matrix. 
% Default: identity matrix.
%
% 'Sigma_nu_delta':
% Covariance for proposal distribution for nu_delta. Early in the burn-in 
% period, this will be replaced with an adaptive covariance matrix. 
% Default: identity matrix.
%
% 'nu_theta_prior_param':
% Hyperparameter c for beta prior distribution on nu_theta, where 
% exp(-exp(nu_theta)) ~ Beta(1,c). Default: 0.4.
%
% 'nu_delta_prior_param':
% Hyperparameter c for beta prior distribution on nu_delta, where 
% exp(-exp(nu_delta)) ~ Beta(1,c). Default: 40.
%
% 'eta':
% The function describing the computer model. Default: Ex_sim.m.
% 
% 'input_cntrl_min':
% Minimum (minima) of control input. Used for normalization. Default: 1.95.
% 
% 'input_cntrl_range':
% Range(s) of control input. Used for normalization. Default: 0.1.
%
% 'input_calib_min':
% Minimum (minima) of calib input. Used for normalization. Default: [0 0].
% Note that the minima should be in the order corresponding to the way the
% input is treated by the model. That is, if the computer model takes
% inputs [t1, t2] and t2 is being treated as state-aware, then in the MCMC
% t2 will be "theta1" and t1 will be "theta2", but still input_calib_min
% should be [t1_min t2_min] (not [theta1_min theta2_min]).
%
% 'input_calib_range':
% Range(s) of calibration input. Used for normalization. Default: [3 6].
% Note that the ranges of the calibration input should be ordered
% corresponding to the way the input is treated by the model -- see the
% discussion of input_calib_min above.
%
% 'output_mean':
% Means of model outputs. Used for standardization. 
% Default: [0.9286    0.7727   21.0030].
% 
% 'output_sd':
% Model standard deviations. Used for standardization.
% Default: [0.0399    0.0776    3.2004]
%
% 'link_fn':
% Link function for theta1. Default: identity.
% 
% 'which_sa'
% Binary array that indicates which calibration inputs are to be treated as
% state-aware. Default: [1 0].
 

%% Parse inputs
p = inputParser;

% Required inputs
p.addRequired('desired_obs',@ismatrix);
p.addRequired('lambda_delta',@isscalar);
p.addRequired('cntrl_input',@ismatrix);

% Varargin
p.addParameter('M',2e4,@isscalar);
p.addParameter('burn_in',1/2,@isscalar);
p.addParameter('nugget',0.05,@isscalar);
p.addParameter('mu_theta',0.5,@ismatrix);
p.addParameter('dim_theta1',1,@isscalar);
p.addParameter('dim_theta2',1,@isscalar);
p.addParameter('lambda_theta_hypers',[0.1 0.1],@ismatrix);
p.addParameter('Sigma_xi',1,@ismatrix);
p.addParameter('Sigma_nu_theta',1,@ismatrix);
p.addParameter('Sigma_nu_delta',1,@ismatrix);
p.addParameter('nu_theta_prior_param',0.4,@isscalar);
p.addParameter('nu_delta_prior_param',40,@isscalar);
p.addParameter('eta',@Ex_sim_compwise);
p.addParameter('input_cntrl_min',1.95,@ismatrix);
p.addParameter('input_cntrl_range',0.1,@ismatrix);
p.addParameter('input_calib_min',[0 0],@ismatrix);
p.addParameter('input_calib_range',[3 6],@ismatrix);
p.addParameter('output_mean',[0.9286    0.7727   21.0030],@ismatrix);
p.addParameter('output_sd',[0.0399    0.0776    3.2004],@ismatrix);
p.addParameter('link_fn',@(x)x,@(x)isa(x,'function_handle'));
p.addParameter('which_sa',[1 0],@ismatrix);

% Parse
p.parse(desired_obs,lambda_delta,cntrl_input,varargin{:});

% Collect inputs
M                    = p.Results.M;
burn_in              = floor(p.Results.burn_in*M);
nugget               = p.Results.nugget;
mu_theta             = p.Results.mu_theta;
dim_theta1           = p.Results.dim_theta1;
dim_theta2           = p.Results.dim_theta2;
lambda_theta_hypers  = p.Results.lambda_theta_hypers;
Sigma_xi             = p.Results.Sigma_xi;
Sigma_nu_theta       = p.Results.Sigma_nu_theta;
Sigma_nu_delta       = p.Results.Sigma_nu_delta;
nu_theta_prior_param = p.Results.nu_theta_prior_param;
nu_delta_prior_param = p.Results.nu_delta_prior_param;
eta                  = p.Results.eta;
input_cntrl_min      = p.Results.input_cntrl_min;
input_cntrl_range    = p.Results.input_cntrl_range;
input_calib_min      = p.Results.input_calib_min;
input_calib_range    = p.Results.input_calib_range;
output_mean          = p.Results.output_mean;
output_sd            = p.Results.output_sd; 
link_fn              = p.Results.link_fn;
which_sa             = p.Results.which_sa;

% Infer useful values
num_out = size(desired_obs,2)  ; % Assume each output is a column 
dim_cntrl = size(cntrl_input,2); % Assume each input is a column

% Create xx, the control inputs including dummy variables
xx = repmat(cntrl_input,num_out,1);
dum_mat = [ eye(num_out - 1) ; zeros(1, num_out -1 ) ] ; % dummy vars
if num_out > 1
    dum_mat = repelem(dum_mat,length(cntrl_input),1);
    xx = [dum_mat xx];
end

% If default Sigma_xi used, make its dim match that of xi
if Sigma_xi == 1
    Sigma_xi = eye(dim_theta2);
end

% If default Sigma_nu_delta used, make its dim correct
if Sigma_nu_delta == 1
    Sigma_nu_delta = eye(dim_cntrl + num_out - 1);
end

% Pack up and leave
settings = struct(...
    'desired_obs',desired_obs,...
    'lambda_delta',lambda_delta,...
    'number_of_iterations',M,...
    'burn_in',burn_in,...
    'cntrl_input',cntrl_input,...
    'cntrl_input_with_dum_vars',xx,...
    'nugget',nugget,...
    'mu_theta',mu_theta,...
    'dim_theta1',dim_theta1,...
    'dim_theta2',dim_theta2,...
    'dim_nu_theta',size(Sigma_nu_theta,1),...
    'dim_output',num_out,...
    'lambda_theta_hypers',lambda_theta_hypers,...
    'Sigma_xi',Sigma_xi,...
    'Sigma_nu_theta',Sigma_nu_theta,...
    'Sigma_nu_delta',Sigma_nu_delta,...
    'nu_theta_prior_param',nu_theta_prior_param,...
    'nu_delta_prior_param',nu_delta_prior_param,...
    'model',eta,...
    'input_cntrl_min',input_cntrl_min,...
    'input_cntrl_range',input_cntrl_range,...
    'input_calib_min',input_calib_min,...
    'input_calib_range',input_calib_range,...
    'output_mean',output_mean,...
    'output_sd',output_sd,...
    'link_fn',link_fn,...
    'which_sa',which_sa);
end
