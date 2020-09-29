function settings = MCMC_dual_calib_settings(...
    sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,...
    des_x,des_y,...
    varargin)
% settings = MCMC_settings(sim_x,sim_t1,sim_t2,sim_y,...
%                          obs_x,obs_t2,obs_y,des_x,des_y,varargin)
% 
% Parameter/Value pairs:
% 'M':
%   Number of samples to draw in the MCMC, includes burn-in. Default: 1e4.
% 'burn_in':
%   Proportion of samples to treat as burn-in. Default: 1/5.
% 'min_x':
%   Vector giving the minimum value of each element of the control input.
%   By default, this is taken to be the minimum of the simulator
%   input provided as sim_x and the observation input obs_x.
% 'range_x':
%   Vector giving the ranges of each element of the control input.
%   By default, this is taken to be the range of the simulator input
%   provided as sim_x and the observation input obs_x.
% 'min_t1':
%   Vector giving the minimum value of each element of the parameter to be
%   calibrated. By default, this is taken to be the minimum of the 
%   simulator input provided as sim_t1.
% 'range_t1':
%   Vector giving the range of each element of calibration parameter.
%   By default, this is taken to be the range of the simulator
%   input provided as sim_t1.
% 'dim_x':
%   Dimensionality of x control input. By default, this is
%   inferred from the x inputs.
% 'dim_t1':
%   Dimensionality of theta1 calibration parameter. By default, this is
%   inferred from the t1 inputs.
% 'dim_t2':
%   Dimensionality of theta2 design parameter. By default, this is
%   inferred from the t2 inputs.
% 'min_t2':
%   Vector giving the minimum value of each element of design parameter.
%   By default, this is taken to be the minimum of the simulator
%   input provided as sim_t2 and the observation input obs_t2.
% 'range_t2':
%   Vector giving the range of each element of design parameter.
%   By default, this is taken to be the range of the simulator
%   input provided as sim_t2 and the observation input obs_t2.
% 'mean_y':
%   Scalar giving the mean value of the output y. By default, this is taken
%   to be the mean of the simulator and observational data sim_y and obs_y.
% 'std_y':
%   Scalar giving the standard deviation of the output y. By default, this
%   is taken to be the standard deviation of the simulator and
%   observational data sim_y and obs_y.
% 'obs_var':
%   Gives the iid gaussian observation error variance. Inputting a
%   scalar value here results in a covariance matrix that is the identity
%   multiplied by that scalar value. Default: scalar 0.05.
% 'des_var':
%   Gives the iid gaussian target outcomes' error variance. Inputting a
%   scalar value here results in a covariance matrix that is the identity
%   multiplied by that scalar value. Default: scalar 0.05.
% 'CTO':
%   Boolean which tells whether the algoritm is being used for CTO. This
%   aids in setting appropriate priors. Default: false.
% 'additional_discrep_mean':
%   Optional additional discrepancy mean. This is primarily intended to
%   be used in the case when CTO is completed after KOH; this additional
%   discrepancy term allows us to include in CTO information about the
%   observation discrepancy learned during KOH. Default 0.
% 'additional_discrep_cov':
%   Optional additional discrepancy covariance. This is primarily intended 
%   to be used in the case when CTO is completed after KOH; this additional
%   discrepancy term allows us to include in CTO information about the
%   observation discrepancy learned during KOH. Default 0.
% 'emulator_use':
%   Determines whether an emulator is used. If not, the emulator mean 
%   function is set to just be the objective function itself, and pairs 
%   this with a 0 covariance fn (by setting marginal precision to Inf).
%   Default: true.
% 'obs_var_est':
%   Boolean value which determines whether a posterior distribution on
%   observation error variance will be explored in MCMC. If not, then the
%   intial value of obs_var is kept constant. Default: false.
% 'obs_var_same':
%   Boolean value which determines whether the observation error variances
%   of the multiple model outputs is constrained to be the same. If set to
%   true, this ensures that the posterior predictive distribution will be
%   centered at the point in the feasible space that is nearest to the
%   target outcome, where all outputs are standardized to have mean 0 and
%   s.d. 1. Default: false.
% 'des_var_est':
%   Boolean value which determines whether a posterior distribution on
%   target error variance will be explored in MCMC. If not, then the
%   intial value of des_var is kept constant. Default: false.
% 'EmulatorMean':
%   Mean function of the emulator. Default: constant 0 (if an emulator is
%   used), example objective function (if emulator==false).
% 'EmulatorCovHypers':
%   Vector of hyperparameters for the power product exponential covariance
%   function of the GP emulator. By default, these are estimated as MLEs.
% 'obs_discrep':
%   Boolean value which tells whether or not to include a discrepancy term
%   for the real observations. Default: true.
% 'des_discrep':
%   Boolean value which tells whether or not to include a discrepancy term
%   for the desired observations. Default: true.
% 'obs_discrep_mean':
%   Function of x and theta2 which gives the prior mean on the observation
%   discrepancy function. Default: constant 0.
% 'obs_rho_lambda':
%   Allows user to specify a value for obs_rho,obs_lambda By default, these
%   are estimated via MCMC.
% 'obs_rho_beta_params':
%   Shape parameters for the beta prior on obs_rho, the characteristic
%   length-scale of the observation discrepancy function. Default: [2,0.4].
% 'obs_lambda_gam_params':
%   Shape/scale parameters for the gamma prior on obs_lambda, the precision
%   of the observation discrepancy function. Default: [5,5].
% 'des_rho_beta_params':
%   Shape parameters for the beta prior on des_rho, the characteristic
%   length-scale of the target discrepancy function. Default: [2,0.4].
% 'des_lambda_gam_params':
%   Shape/scale parameters for the gamma prior on des_lambda, the precision
%   of the target discrepancy function. Default: [150,4].
% 'obs_final_size':
%   Final size of the set of observations, if sequential DoE is to be used.
%   By default this is taken to be the size(obs_x,1), so that no new
%   observations are created in the course of DCTO.
% 'true_phenomenon':
%   If new observations are made during the course of DCTO (as in
%   sequential DoE), then the high fidelity model must be provided here,
%   from which the new observations will be obtained. If sequential DoE is
%   used and no function is provided here, an exception will occur.
% 'true_obs_var':
%   If new observations are made during the course of DCTO (as in
%   sequential DoE), then the true observation error variance must be
%   provided here, so that it can be applied to new observations.
% 'obs_discrep_use_MLEs':
%   Boolean value which tells whether or not to estimate the observation
%   discrepancy covariance hyperparameters via maximum likelihood
%   estimation periodically during the burn-in process using a point
%   estimate of the calibration parameter. Default: false.
% 'modular':
%   Boolean value which tells whether or not to use modularized version of
%   the model, to protect the traditional calibration from being influenced
%   by the desired observations. Default: false.
% 'doplot':
%   true       - (Default) Update scatterplot every ten draws.
%   false      - No plots during MCMC.
% 'verbose':
%   true       - (Default) Output verbose progress descriptions.
%   false      - Output minimal progress descriptions.


%% Parse inputs
p = inputParser;
p.addRequired('sim_x',@ismatrix);
p.addRequired('sim_t1',@ismatrix);
p.addRequired('sim_t2',@ismatrix);
p.addRequired('sim_y',@ismatrix);
p.addRequired('obs_x',@ismatrix);
p.addRequired('obs_t2',@ismatrix);
p.addRequired('obs_y',@ismatrix);
p.addRequired('des_x',@ismatrix);
p.addRequired('des_y',@ismatrix);
p.addParameter('M',1e4,@isscalar);
p.addParameter('burn_in',1/5,@(x) x>0 && x<1);
p.addParameter('min_x','Default',@isnumeric);
p.addParameter('range_x','Default',@isnumeric);
p.addParameter('min_t1','Default',@isnumeric);
p.addParameter('range_t1','Default',@isnumeric);
p.addParameter('dim_x','Default',@isscalar);
p.addParameter('dim_t1','Default',@isscalar);
p.addParameter('dim_t2','Default',@isscalar);
p.addParameter('min_t2','Default',@isnumeric);
p.addParameter('range_t2','Default',@isnumeric);
p.addParameter('mean_y','Default',@isnumeric);
p.addParameter('std_y','Default',@isnumeric);
p.addParameter('obs_var',0.05,@ismatrix);
p.addParameter('des_var',0,@ismatrix);
p.addParameter('CTO',false,@islogical);
p.addParameter('additional_discrep_mean',@(x)0,...
    @(h)isa(h,'function_handle'));
p.addParameter('additional_discrep_cov','Default',...
    @(h)isa(h,'function_handle'));
p.addParameter('emulator_use',true,@islogical);
p.addParameter('obs_var_est',false,@islogical);
p.addParameter('obs_var_same',false,@islogical);
p.addParameter('des_var_est',false,@islogical);
p.addParameter('EmulatorMean','Default',@(h)isa(h,'function_handle'));
p.addParameter('EmulatorCovHypers','Default',@ismatrix);
p.addParameter('obs_discrep',true,@islogical);
p.addParameter('des_discrep',true,@islogical);
p.addParameter('obs_discrep_mean','Default',@(h)isa(h,'function_handle'));
p.addParameter('obs_rho_lambda','Default',@isnumeric);
p.addParameter('obs_rho_beta_params','Default',@isnumeric);
p.addParameter('obs_lambda_gam_params','Default',@isnumeric);
p.addParameter('des_rho_beta_params','Default',@isnumeric);
p.addParameter('des_lambda_gam_params','Default',@isnumeric);
p.addParameter('obs_final_size',size(obs_x,1),@isscalar);
p.addParameter('true_phenomenon',...
    @(varargin)error('Error: true phenomenon not supplied'),...
    @(h)isa(h,'function_handle'));
p.addParameter('true_obs_var',0.05,@isnumeric);
p.addParameter('obs_discrep_use_MLEs',false,@islogical);
p.addParameter('modular',false,@islogical);
p.addParameter('doplot',true,@islogical);
p.addParameter('verbose',true,@islogical);
p.parse(sim_x,sim_t1,sim_t2,sim_y,obs_x,obs_t2,obs_y,des_x,des_y,...
    varargin{:});

%% Collect inputs
M = p.Results.M;
burn_in = floor(p.Results.burn_in*M);
min_x = p.Results.min_x;
range_x = p.Results.range_x;
min_t1 = p.Results.min_t1;
range_t1 = p.Results.range_t1;
dim_x = p.Results.dim_x; 
dim_t1 = p.Results.dim_t1; 
dim_t2 = p.Results.dim_t2;
min_t2 = p.Results.min_t2;
range_t2 = p.Results.range_t2;
mean_y = p.Results.mean_y;
std_y = p.Results.std_y;
obs_var = p.Results.obs_var;
des_var = p.Results.des_var;
CTO = p.Results.CTO;
additional_discrep_cov = p.Results.additional_discrep_cov;
additional_discrep_mean = p.Results.additional_discrep_mean;
emulator_use = p.Results.emulator_use;
obs_var_est = p.Results.obs_var_est;
obs_var_same = p.Results.obs_var_same;
des_var_est = p.Results.des_var_est;
EmulatorMean = p.Results.EmulatorMean;
EmulatorCovHypers = p.Results.EmulatorCovHypers;
obs_discrep = p.Results.obs_discrep;
des_discrep = p.Results.des_discrep;
obs_discrep_mean = p.Results.obs_discrep_mean;
obs_rho_lambda = p.Results.obs_rho_lambda;
obs_rho_beta_params = p.Results.obs_rho_beta_params;
obs_lambda_gam_params = p.Results.obs_lambda_gam_params;
des_rho_beta_params = p.Results.des_rho_beta_params;
des_lambda_gam_params = p.Results.des_lambda_gam_params;
obs_final_size = p.Results.obs_final_size;
true_phenomenon = p.Results.true_phenomenon;
true_obs_var = p.Results.true_obs_var;
obs_discrep_use_MLEs = p.Results.obs_discrep_use_MLEs;
modular = p.Results.modular;
doplot = p.Results.doplot;
verbose = p.Results.verbose;

%% Todo
if ~(isequal(obs_rho_lambda,'Default'))
    error('Error: obs_rho_lambda input parameter not implemented yet!');
end


%% Some useful definitions
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x) ./ (1+exp(x));

%% If emulator not used, then give empty simulator arrays appropriate dims
% Avoids headaches with concatenation.
if ~emulator_use
    if numel(sim_x)+numel(sim_t1)+numel(sim_t2) > 0
        error(...
            'Error: emulator_use set to false, but simulations supplied.');
    end
    sim_x = zeros(0,size(min_x,2));
    sim_t1 = zeros(0,size(min_t1,2));
    sim_t2 = zeros(0,size(min_t2,2));
    sim_y = zeros(0,size(obs_y,2));
end

%% If targets not used, then give empty target arrays appropriate dims
% Avoids headaches with concatenation
if numel(des_x) == numel(des_y) && numel(des_x) == 0
    des_x = zeros(0,size(min_x,2));
    des_y = zeros(0,size(obs_y,2));
end


%% Normalize inputs
x = [sim_x ; obs_x] ;
if isequal(min_x,'Default'), min_x = min(x) ; end
if isequal(range_x,'Default'), range_x = range(x) ; end
if sum(size(sim_x))>0, sim_x_01 = (sim_x - min_x) ./ range_x ;
else, sim_x_01 = sim_x ; end
obs_x_01 = (obs_x - min_x) ./ range_x ;
if sum(size(des_x))>0, des_x_01 = (des_x - min_x) ./ range_x ;
else, des_x_01 = des_x ; end

if isequal(min_t1,'Default'), min_t1 = min(sim_t1) ; end
if isequal(range_t1,'Default'), range_t1 = range(sim_t1) ; end
if numel(sim_t1)>0, sim_t1_01 = (sim_t1 - min_t1) ./ range_t1 ;
else, sim_t1_01 = sim_t1; end

t2 = [sim_t2 ; obs_t2 ] ;
if isequal(min_t2,'Default'), min_t2 = min(t2) ; end
if isequal(range_t2,'Default'), range_t2 = range(t2) ; end
sim_t2_01 = (sim_t2 - min_t2) ./ range_t2 ; 
obs_t2_01 = (obs_t2 - min_t2) ./ range_t2 ;


%% Infer some useful values
if isequal(dim_x,'Default')
    dim_x = numel(min_x);
end
if isequal(dim_t1,'Default')
    dim_t1 = numel(min_t1);
end
if isequal(dim_t2,'Default')
    dim_t2 = numel(min_t2);
end
dim_y = size(obs_y,2);


%% Standardize outputs
if isequal(mean_y,'Default'), mean_y = mean(sim_y) ; end
if isequal(std_y,'Default'), std_y = std(sim_y) ; end
if numel(sim_y) == 0, sim_y = repmat(obs_y,0,1) ; end
sim_y_std = (sim_y - mean_y) ./ std_y ; 
obs_y_std = (obs_y - mean_y) ./ std_y ;
des_y_std = (des_y - mean_y) ./ std_y ; 



%% Set emulator mean and covariance hyperparameters
% Determines whether an emulator is used. If not, the emulator 
% mean function is set to just be the objective
% function itself, and pairs this with a 0
% covariance function (by setting marginal precision to Inf).
if emulator_use
    if isequal(EmulatorMean,'Default')
        mean_sim = @(a) zeros(size(a,1),dim_y); % Emulator mean
    else
        mean_sim = EmulatorMean;
    end
    if isequal(EmulatorCovHypers,'Default')
        EmulatorCovHypers = zeros(dim_x+dim_t1+dim_t2+1,dim_y);
        mean_sim_vals = mean_sim([sim_x_01,sim_t1_01,sim_t2_01]);
        % Loop through the outputs of the function
        for ii = 1:dim_y
            % Define function for minimization
            f = @(rl) ...
                -logmvnpdf(sim_y_std(:,ii)',...
                mean_sim_vals(:,ii)',...
                gp_cov(rl(1:(end-1)),...
                [sim_x_01 sim_t1_01 sim_t2_01],...
                [sim_x_01 sim_t1_01 sim_t2_01],...
                rl(end),false) + ...
                1e-4*eye(size(sim_t1_01,1)));
            % Perform minimization
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = [zeros(1,size([sim_x sim_t1 sim_t2],2)) 0];
            ub = [ones(1,size([sim_x sim_t1 sim_t2],2)) Inf];
            x0 = [.5*ones(1,size([sim_x sim_t1 sim_t2],2)) 1];
            [inp,~,~,~] = fmincon(f,x0,A,b,Aeq,beq,lb,ub);
            EmulatorCovHypers(:,ii) = inp ;
        end
    end
else
    if isequal(EmulatorMean,'Default')
        mean_sim = @(a) dual_calib_example_fn(...
            a(:,1:length(min_x)),min_x,range_x,a(:,((length(min_x)+1):(length(min_x)+length(min_t1)))),min_t1,range_t1,a(:,( (length(min_x) + length(min_t1) + 1 ):end)),min_t2,range_t2,...
            mean_y,std_y,0,true);
    else
        mean_sim = EmulatorMean;
    end
    EmulatorCovHypers = ...
        repmat([.5* ones(1,dim_x+dim_t1+dim_t2) Inf]',1,dim_y);
end

%% Set observation discrepancy mean
if isequal(obs_discrep_mean,'Default')
    mean_obs = @(x,t) zeros(size(x,1),dim_y);
else
    mean_obs = obs_discrep_mean;
end

%% Set theta priors
log_theta1_prior = @(t) 0 ; % Uniform prior
log_theta2_prior = @(t) 0 ; % Uniform prior


%% Set theta proposal distributions and initial values
% We'll draw logit-transformed theta from a normal centered at the
% logit-transformed previous draw.
theta1_proposal = @(t,S) logit_inv(mvnrnd(logit(t),S)); 
% theta1_proposal = @(t,S) rand(1,length(t)); 
theta2_proposal = @(t,S) logit_inv(mvnrnd(logit(t),S)); 
% Since this proposal is not symmetric, we need to use full
% Metropolis-Hastings rather than just Metropolis. So here is the log MH
% correction for the lack of symmetry.
theta1_prop_log_mh_correction = ...
    @(t_s,t) sum(log(t_s)+log(1-t_s)-log(t)-log(1-t)); 
% theta1_prop_log_mh_correction = ...
%     @(t_s,t) 0;
theta2_prop_log_mh_correction = ...
    @(t_s,t) sum(log(t_s)+log(1-t_s)-log(t)-log(1-t));
% Set initial values and initial covariance matrices for proposals
theta1_init = rand(dim_t1,1);
theta2_init = rand(dim_t2,1);
Sigma_theta1 = eye(size(theta1_init,1));
Sigma_theta2 = eye(size(theta2_init,1));


%% Set real and target discrepancy rho and lambda prior distributions
if isequal(obs_rho_beta_params,"Default")
    obs_rho_beta_params = [2,0.4];
end
if isequal(obs_lambda_gam_params,"Default")
    obs_lambda_gam_params = [5,5];
end
if isequal(des_rho_beta_params,"Default")
    des_rho_beta_params = [2,0.4];
end
if isequal(des_lambda_gam_params,"Default")
    des_lambda_gam_params = [150,4];
end
log_obs_rho_prior  = @(r) sum(log( betapdf(...
    r,obs_rho_beta_params(1),obs_rho_beta_params(2)) ));
log_obs_lambda_prior = @(ld) log( gampdf(...
    ld,obs_lambda_gam_params(1),obs_lambda_gam_params(2)) );
log_des_rho_prior  = @(r) sum(log( betapdf(...
    r,des_rho_beta_params(1),des_rho_beta_params(2)) ));
log_des_lambda_prior = @(ld) log( gampdf(...
    ld,des_lambda_gam_params(1),des_lambda_gam_params(2)) );

%% Set real and target error variance prior distributions
if CTO
%     log_obs_var_prior_fn = @(sig2) log( gampdf(...
%         sig2,4,.125) );
    log_obs_var_prior_fn = @(sig2) log( exppdf(...
        sig2,.001) );
    log_des_var_prior_fn = @(sig2) error('Error: Is KOH, CTO or DCTO?');
else
%     log_obs_var_prior_fn = @(sig2) -log(sig2);
%     log_des_var_prior_fn = @(sig2) -8 * log( sig2 );
    log_obs_var_prior_fn = @(sig2) log( gampdf(...
    sig2,1,1) );
    log_des_var_prior_fn = @(sig2) log( gampdf(...
    sig2,4,.125) );
end



%% Set real and desired discrepancy rho and lambda proposal distributions
% if isequal(obs_rho_lambda,'Default')
% We'll draw logit-transformed rho, and lambda from normals centered at
% log-transformed previous draw.
rho_proposal = @(r,S) logit_inv(mvnrnd(logit(r),S)); 
lambda_proposal = @(lam,S) exp(mvnrnd(log(lam),S)); 
% Since this proposal is not symmetric, we need to use full
% Metropolis-Hastings rather than just Metropolis. So here is the log MH
% correction for the lack of symmetry.
rho_prop_log_mh_correction = ...
    @(r_s,r) sum(log(r_s)+log(1-r_s)-log(r)-log(1-r));
lambda_prop_log_mh_correction = ...
    @(lam_s,lam) sum(log(lam_s)-log(lam));
% Set initial values and initial covariance matrices for proposals
if obs_discrep
    obs_rho_init = betarnd(...
        obs_rho_beta_params(1),obs_rho_beta_params(2),dim_x+dim_t2,dim_y);
    obs_lambda_init = gamrnd(obs_lambda_gam_params(1),...
        obs_lambda_gam_params(2),1,dim_y);
else
    obs_rho_init = .5*ones(dim_x+dim_t2,dim_y);
    obs_lambda_init = Inf * ones(1,dim_y);
end
if des_discrep
    des_rho_init = betarnd(...
        des_rho_beta_params(1),des_rho_beta_params(2),dim_x,dim_y);
    des_lambda_init = gamrnd(des_lambda_gam_params(1),...
        des_lambda_gam_params(2),1,dim_y);
else
    des_rho_init = .5*ones(size(des_x,2),dim_y);
    des_lambda_init = Inf * ones(1,dim_y);
end

% Get initial proposal variances
obs_Sigma_rho = repmat(eye(size(obs_rho_init,1)),1,1,dim_y);
obs_Sigma_lambda = ones(1,dim_y);
des_Sigma_rho = repmat(eye(size(des_rho_init,1)),1,1,dim_y);
des_Sigma_lambda = ones(1,dim_y);

%% Set up observation and target err var arrays, proposals, log mh correctn
% We need an array of obs and des variances, one for each model output
if isscalar(obs_var)
    obs_var = ones(1,dim_y) * obs_var ; 
end
if isscalar(des_var)
    des_var = ones(1,dim_y) * des_var ; 
end

% We need proposal distributions for obs_var and des_var: lognormal
% distribution on the log-transform of obs_var
obs_var_proposal = @(sig,S) exp(randn * sqrt(S) + log(sig)); 
des_var_proposal = @(sig,S) exp(randn * sqrt(S) + log(sig)); 
% Now we need the log metropolis-hastings corrections for these proposals
obs_var_prop_log_mh_correction = ...
    @(sig_s,sig) sum(log(sig_s)-log(sig));
des_var_prop_log_mh_correction = ...
    @(sig_s,sig) sum(log(sig_s)-log(sig));

% Get initial proposal variances
obs_var_Sigma = ones(1,dim_y) ; 
des_var_Sigma = ones(1,dim_y) ; 

%% Make additional covariance matrix (used for adding discrepancy cov)
if isequal(additional_discrep_cov,'Default')
    additional_discrep_cov = @(x,theta) ...
        zeros(max(1,size(x,1)),max(1,size(x,1)),dim_y) ; 
end

%% Make obs_var_est logical the appropriate dim
if isscalar(obs_var_est) && dim_y ~= 1
    obs_var_est = repmat(obs_var_est,1,dim_y);
end

%% Pack up the settings structure
settings = struct(...
    'M',M,...
    'burn_in',burn_in,...
    'sim_x',sim_x_01,...
    'sim_t1',sim_t1_01,...
    'sim_t2',sim_t2_01,...
    'sim_y',sim_y_std,...
    'obs_x',obs_x_01,...
    'obs_t2',obs_t2_01,...
    'obs_y',obs_y_std,...
    'des_x',des_x_01,...
    'des_y',des_y_std,...
    'min_x',min_x,...
    'range_x',range_x,...
    'min_t1',min_t1,...
    'range_t1',range_t1,...
    'min_t2',min_t2,...
    'range_t2',range_t2,...
    'mean_y',mean_y,...
    'std_y',std_y,...
    'obs_var',obs_var,...
    'des_var',des_var,...
    'additional_discrep_cov',additional_discrep_cov,...
    'additional_discrep_mean',additional_discrep_mean,...
    'mean_sim',mean_sim,...
    'emulator_use',emulator_use,...
    'emulator_rho',EmulatorCovHypers(1:end-1,:),...
    'emulator_lambda',EmulatorCovHypers(end,:),...
    'mean_obs',mean_obs,...
    'theta1_proposal',theta1_proposal,...
    'theta2_proposal',theta2_proposal,...
    'rho_proposal',rho_proposal,...
    'lambda_proposal',lambda_proposal,...
    'obs_var_est',obs_var_est,...
    'obs_var_same',obs_var_same,...
    'des_var_est',des_var_est,...
    'obs_var_proposal',obs_var_proposal,...
    'des_var_proposal',des_var_proposal,...
    'theta1_prop_log_mh_correction',theta1_prop_log_mh_correction,...
    'theta2_prop_log_mh_correction',theta2_prop_log_mh_correction,...
    'rho_prop_log_mh_correction',rho_prop_log_mh_correction,...
    'lambda_prop_log_mh_correction',lambda_prop_log_mh_correction,...
    'obs_var_prop_log_mh_correction',obs_var_prop_log_mh_correction,...
    'des_var_prop_log_mh_correction',des_var_prop_log_mh_correction,...
    'des_rho_init',des_rho_init,...
    'des_lambda_init',des_lambda_init,...
    'des_Sigma_rho',des_Sigma_rho,...
    'des_Sigma_lambda',des_Sigma_lambda,...
    'obs_var_Sigma',obs_var_Sigma,...
    'des_var_Sigma',des_var_Sigma,...
    'obs_rho_init',obs_rho_init,...
    'obs_lambda_init',obs_lambda_init,...
    'obs_Sigma_rho',obs_Sigma_rho,...
    'obs_Sigma_lambda',obs_Sigma_lambda,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'log_obs_rho_prior',log_obs_rho_prior,...
    'log_obs_lambda_prior',log_obs_lambda_prior,...
    'log_des_rho_prior',log_des_rho_prior,...
    'log_des_lambda_prior',log_des_lambda_prior,...
    'log_obs_var_prior_fn',log_obs_var_prior_fn,...
    'log_des_var_prior_fn',log_des_var_prior_fn,...
    'theta1_init',theta1_init,...
    'theta2_init',theta2_init,...
    'theta1_prop_cov',Sigma_theta1,...
    'theta2_prop_cov',Sigma_theta2,...
    'log_theta1_prior',log_theta1_prior,...
    'log_theta2_prior',log_theta2_prior,...
    'obs_discrep',obs_discrep,...
    'des_discrep',des_discrep,...
    'obs_final_size',obs_final_size,...
    'true_phenomenon',true_phenomenon,...
    'true_obs_var',true_obs_var,...
    'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
    'modular',modular,...
    'doplot',doplot,...
    'verbose',verbose);

end