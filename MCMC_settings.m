function settings = MCMC_settings (M,desired_obs,which_outputs)


%% MCMC settings
burn_in = ceil(M/5) ; % burn-in
% Covariance parameter settings, found by optimization routine:
Rho_lam_optimum  = [0.655344235568109   0.931941001705886 ...
    0.960653924901867   0.991953924787049   0.017385994893994];
omega  = Rho_lam_optimum(1:2);
rho    = Rho_lam_optimum(3:4);
lambda = Rho_lam_optimum(5);
num_out = length(desired_obs);

%% Proposal density
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x)./(1+exp(x));
prop_density = @(x,Sigma) logit_inv(mvnrnd(logit(x),Sigma)); % Normal
Sigma = [.5 0 ; 0 .5]; % Initial variance for prop_density
% Upper and lower bounds for theta (if applicable)
LB = [ .1 10 ] ; UB = [.6 25] ;
% Nugget size for computational stability of covariance matrices
nugsize = @(Covmat) 1e-4 ; % Why make it a function? So later it can be 
                           % made fancier.
                           
%% Set prior for sigma^2_y: 1/sigma^2_y, and log version
sigma2 = 1; % initial value
sigma2_prior = @(sigma2) 1/prod(sigma2);
log_sigma2_prior = @(sigma2) -log(prod(sigma2));
sigma2_prop_density = @(x,s) normrnd(x,s);
Sigma_sig = 2 ;
% Piecewise version:
sigma2 = rand(1,num_out)*10;
sigma2_prop_density = @(x,s) normrnd(x,s);
Sigma_sig = 4 * ones(1,num_out);
% Joint multivariate draw version:
sigma2=rand(1,num_out)*20;
sigma2_prop_density = @(x,s) exp(mvnrnd(log(x),s));
log_sig_mh_correction = @(sig_s,sig) log(prod(sig_s)) - log(prod(sig));
Sigma_sig = eye(num_out);

%% Set prior for theta
Cost_lambda = 100;
Cost = @(t,Cost_lambda) Cost_lambda * norm(t)^2;
theta_prior = @(theta) exp(-Cost(theta,Cost_lambda));
log_theta_prior = @(theta) -Cost(theta,Cost_lambda);


%% Package proposal density
proposal.density = prop_density; 
proposal.Sigma = Sigma;
proposal.sigma2_prop_density = sigma2_prop_density;
proposal.Sigma_sig = Sigma_sig;

%% MH correction for using log-normal proposal
log_mh_correction = @(theta_s,theta) log(prod(theta_s)*prod(1-theta_s))-...
    log(prod(theta)*prod(1-theta));


%% Load data and get initial theta value
fprintf('Reading data from .xlsx...\n')
raw_dat = xlsread('fe_results.xlsx');
indx = [1 2 3]; % This tells us which columns of raw_dat we need.
for ii = 1:length(which_outputs) % This loop will set indx appropriately.
    if which_outputs(ii) indx = [indx 3+ii] ; end
end
raw_dat = raw_dat(:,indx);
% raw_dat is assumed to include one observation per row, with input columns
% preceding output columns (with no headers).
num_calib = 2 ; % This is the number of calibration parameters in the data,
                % which are assumed to be the columns immediately preceding
                % the output columns in raw_dat.
num_cntrl = size(raw_dat,2) - num_out - num_calib + (num_out>1) ; 
% num_cntrl is the number of control inputs in raw_dat, plus one if the
% output is multivariate; the extra one is for a dummy input.
init_theta = rand(1,2) 

%% Rescale inputs, standardize outputs
tdat = Tdat(raw_dat,num_out); % rescaling inputs, standardizing outputs. 
                              % Tdat also adds a dummy control variable
                              % with num_out levels.
fprintf('done.\n\n')
sim_xt = tdat.input ;
eta = tdat.output;
% Input mins, input ranges, output means, and output sd's are stored so
% that the rescaling and standardization can be reversed later.
sim_cntrl_input_mins = tdat.input_mins(1:num_cntrl) ;
sim_cntrl_input_ranges = tdat.input_ranges(1:num_cntrl) ;
sim_calib_input_mins = tdat.input_mins(num_cntrl+1:end) ;
sim_calib_input_ranges = tdat.input_ranges(num_cntrl+1:end) ;
sim_output_means = mean(tdat.output_means') ;
% Note that in the above line, the mean is taken across temperatures. This
% is because the standardization was carried out by individual temperature
% setting, whereas since our desired data will be constant across
% temperature, we want a single value for all temperature settings. This
% fudge should be harmless, in that the values do not differ much across
% temperature. Sd's will be treated similarly in the following line.
sim_output_sds = mean(tdat.output_sds') ;
% Assuming uniform prior for theta: find LB,UB and rescale them
LB = (LB - sim_calib_input_mins)./sim_calib_input_ranges ;
UB = (UB - sim_calib_input_mins)./sim_calib_input_ranges ;

%% Prepare field observations, control settings, and variance for MCMC
% First, we must standardize them to be on the same scale as the simulation
% observations.
y = (desired_obs - sim_output_means) ./ sim_output_sds ;
% Now, pair these standardized observations with appropriate control
% settings. Notice that we are here assuming only one field observation,
% which is constant across control settings (other than dummy variable).
obs_x = unique(tdat.input(:,1:num_cntrl),'rows');
y = repelem(y,size(obs_x,1)/length(y))' ;


%% Set uniform prior for theta
out_of_range = @(theta) theta < LB | theta > UB ; 

settings = struct('M',M,'burn_in',burn_in,'sim_xt',sim_xt,...
    'eta',eta,'obs_x',obs_x,'y',y,'sigma2',sigma2,...
    'log_sigma2_prior',log_sigma2_prior,'out_of_range',out_of_range,...
    'init_theta',init_theta,'omega',omega,'rho',rho,'lambda',lambda,...
    'proposal',proposal,'nugsize',nugsize,'num_out',num_out,...
    'log_sig_mh_correction',log_sig_mh_correction,...
    'log_mh_correction',log_mh_correction,...
    'output_sds',tdat.output_sds,'output_means',tdat.output_means,...
    'log_theta_prior',log_theta_prior,...
    'Cost_lambda',Cost_lambda,...
    'which_outputs',which_outputs);

end