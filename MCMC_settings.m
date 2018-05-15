function settings = MCMC_settings (M,desired_obs,sim_x,sim_t,sim_y,...
    which_outputs, Rho_lam_optimum)

%% Infer useful values
num_cal = size(sim_t,2);
num_out = size(sim_y,2);
num_cntrl = size(sim_x,2) + num_out - 1 ; % The num_out - 1 is for dum vars
num_cntrl_wod = size(sim_x,2) ; % Without dummy vars

%% MCMC settings
burn_in = ceil(M/5) ; % burn-in
% Covariance parameter settings, found by optimization routine:
if Rho_lam_optimum == 0
    fprintf('No rho,lambda values specified; commencing ML estimation.\n');
    Rho_lam_optimum = opt_rho_lambda(sim_xt,sim_dat,num_cal,...
        rand(1,num_cntrl),rand(1,num_cal),gamrnd(5,5));
end
% Calculated optimum for wind turbine application:
%Rho_lam_optimum  = [   0.935753521438069   0.650946653103927...
%    0.673593619101900  0.479684392594821   0.967330479380613...
%    0.015203646313917 ] ;
% Calculated optimum for example simulation:
%Rho_lam_optimum  = [   0.280981573480363   0.999189406633873...
%    0.600440750045477  0.719652153362981   0.102809702497319...
%    0.000837772517865 ] ;
omega  = Rho_lam_optimum(1:num_cntrl);
rho    = Rho_lam_optimum(num_cntrl+1:num_cntrl+num_cal);
lambda = Rho_lam_optimum(num_cntrl+num_cal+1);
%num_out = length(desired_obs);
% Need different omega in case not all three outputs used
if sum(which_outputs) ~= 3
    error('To use <3 outputs, sort out manually what to do with omega');
    %omega = omega(which_outputs==1);
end

%% Proposal density
logit = @(x) log(x./(1-x));
logit_inv = @(x) exp(x)./(1+exp(x));
prop_density = @(x,Sigma) logit_inv(mvnrnd(logit(x),Sigma)); % Normal
Sigma = [.5 0 ; 0 .5]; % Initial variance for prop_density
% Upper and lower bounds for theta (if applicable)
LB = min(sim_t) ; UB = max(sim_t) ;
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
Cost_lambda = 0;
Cost = @(t,Cost_lambda) Cost_lambda * norm(t)^2;
theta_prior = @(theta,Cost_lambda) exp(-Cost(theta,Cost_lambda));
log_theta_prior = @(theta,Cost_lambda) -Cost(theta,Cost_lambda);


%% Package proposal density
proposal.density             = prop_density; 
proposal.Sigma               = Sigma;
proposal.sigma2_prop_density = sigma2_prop_density;
proposal.Sigma_sig           = Sigma_sig;

%% MH correction for using log-normal proposal
log_mh_correction = @(theta_s,theta) log(prod(theta_s)*prod(1-theta_s))-...
    log(prod(theta)*prod(1-theta));


%% Load data and get initial theta value
% fprintf('Reading data from .xlsx...\n')
% raw_dat = xlsread('fe_results.xlsx');
raw_dat = [sim_x sim_t sim_y];
indx = 1:num_out; % This will help tell which columns of raw_dat we need.
for ii = 1:length(which_outputs) % This loop will set indx appropriately.
    if which_outputs(ii) indx = [indx num_cal+num_cntrl_wod+ii] ; end
end
raw_dat = raw_dat(:,indx);

% Get initial theta val
init_theta = rand(1,num_cal) 

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
obs_x = unique(tdat.input(:,1:num_cntrl),'rows','stable');
y = repelem(y,size(obs_x,1)/length(y))' ;


%% Set uniform prior for theta
out_of_range = @(theta) theta < LB | theta > UB ; 

settings = struct(...
    'M',M,...
    'burn_in',burn_in,...
    'sim_xt',sim_xt,...
    'num_cntrl',num_cntrl,...
    'num_cntrl_wod',num_cntrl_wod,...
    'num_cal',num_cal,...
    'num_out',num_out,...
    'eta',eta,...
    'obs_x',obs_x,...
    'y',y,...
    'sigma2',sigma2,...
    'log_sigma2_prior',log_sigma2_prior,...
    'out_of_range',out_of_range,...
    'init_theta',init_theta,...
    'omega',omega,...
    'rho',rho,...
    'lambda',lambda,...
    'proposal',proposal,...
    'nugsize',nugsize,...
    'log_sig_mh_correction',log_sig_mh_correction,...
    'log_mh_correction',log_mh_correction,...
    'input_cntrl_mins',min(sim_x),...
    'input_calib_mins',min(sim_t),...
    'input_cntrl_ranges',range(sim_x),...
    'input_calib_ranges',range(sim_t),...
    'output_sds',tdat.output_sds,...
    'output_means',tdat.output_means,...
    'log_theta_prior',log_theta_prior,...
    'Cost_lambda',Cost_lambda,...
    'which_outputs',which_outputs,...
    'desired_obs',desired_obs);

end