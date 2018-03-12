% Workspace: Optimize rho, omega for average lambda

clc; clear all; close all;

%% Get data
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');

fprintf('Reading data from .xlsx...\n')
raw_dat = xlsread("NSF DEMS\Phase 1\stored_data\fe_results.xlsx");
tdat = Tdat(raw_dat,3); % rescaling inputs, standardizing outputs
fprintf('done.\n\n')

%% Covariance parameter settings
% Found by optimization routine below:
Rho_lam_optimum  = [0.655344235568109   0.931941001705886 ...
    0.960653924901867   0.991953924787049   0.017385994893994];
omega  = Rho_lam_optimum(1:2);
rho    = Rho_lam_optimum(3:4);
lambda = Rho_lam_optimum(5);

%% Narrow the data set for convenience
raw_dat_s = raw_dat(mod(1:size(raw_dat,1),3)==0,:);
tdat = Tdat(raw_dat_s,3);

%% Assign variables
theta = [.5 .5];
obs_x = [];
obs_y = [];
sim_x = tdat.input(:,1:2);
sim_t = tdat.input(:,3:4);
sim_y = tdat.output;
cntrl_input = sim_x;
calib_input = sim_t;
output = sim_y;
c     = 10^3.625;

%% Generate some observed data
defl_sd = 75/4;
defl_mean = 75/2;
rot_sd = .09/4;
rot_mean = .09/2;
cost_sd = 264/4;
cost_mean = 264/2;
num_obs = 5; % This is how many observations we want
defl_vals = randn(num_obs,1)*defl_sd + defl_mean;
rot_vals = randn(num_obs,1)*rot_sd + rot_mean;
cost_vals = randn(num_obs,1)*cost_sd + cost_mean;
% Build observation input matrix
obs_x = repmat(unique(raw_dat(:,1)),num_obs,1);
ones_vec = ones(size(obs_x));
obs_x = [0*ones_vec obs_x ; .5 * ones_vec obs_x ; ones_vec obs_x ];
% Build the response vector
num_temps = length(unique(raw_dat(:,1)));
obs_y = [ repelem(defl_vals,num_temps,1) ; ...
    repelem(rot_vals,num_temps,1) ; repelem(cost_vals,num_temps,1)];


%% Run the likelihood function
likelihood = Lorl(omega,rho,lambda,obs_x,rho,obs_y,...
    tdat.input(:,1:2),tdat.input(:,3:4),tdat.output,c)

%% Run the log likelihood function
log_likelihood = logLorl(omega,rho,lambda,obs_x,rho,obs_y,...
    tdat.input(:,1:2),tdat.input(:,3:4),tdat.output)


%% Set prior distribution of lambda
a = 5;
b = 1/5;
pi_lambda = @(L) gampdf(L,a,b);
% Support is positive reals

%% Set prior distribution of omega,rho
pi_Rho = @(Rho) prod((ones(4,1) - Rho).^(-.9));
% Support is [0,1]

%% Define integrand
intgrnd = @(Rho,L,obs_x,theta,obs_y,cntrl_input,calib_input,output,c) ...
    Lorl(Rho(1:2)',Rho(3:4)',L,obs_x,theta,obs_y,cntrl_input, ...
    calib_input,output,c) .* pi_lambda(L);

%% Define integral (as function of omega, rho)
intgrl = @(Rho,obs_x,theta,obs_y,cntrl_input,calib_input,output,c) ...
    integral(@(L) intgrnd(Rho,L,obs_x,theta,obs_y,cntrl_input, ...
    calib_input,output,c),0.01,2.25);

%% Define optimization fn of Rho over the integral (given data, theta, etc)
options = optimset('Display','iter','PlotFcns',@optimplotfval,...
    'TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',500);
Rho_opt = @(obs_x,theta,obs_y,cntrl_input,...
    calib_input,output,c,options)...
    fminsearch( @(Rho) -intgrl(Rho,obs_x,theta,obs_y,cntrl_input, ...
    calib_input,output,c),[omega' ; rho' ],options);

%% Perform optimization
[Rho_optimum,fval,exitflag,outpt] = Rho_opt(obs_x,theta,obs_y,...
    cntrl_input,calib_input,output,c,options)

%intgrnd(Rho,L,obs_x,theta,obs_y,cntrl_input,calib_input,output)
%Lorl(Rho(1:2),Rho(3:4),L,obs_x,theta,obs_y,cntrl_input,calib_input,output)

%% Now let's try without integrating lambda out.
options = optimset('Display','iter','PlotFcns',@optimplotfval,...
    'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',500);
%c=10^3.625;
Rho_lam_init = [omega rho lambda]; %Set initial value
lb = [ 0 0 0 0 0 ];
ub = [ 1 1 1 1 Inf ] ;
Rho_lambda_opt = @(obs_x,theta,obs_y,cntrl_input,...
    calib_input,output,options)...
    fmincon(@(Rho_lam) -logLorl(Rho_lam(1:2),Rho_lam(3:4),...
    Rho_lam(5),obs_x,theta,obs_y,cntrl_input,calib_input,output),...
    Rho_lam_init,[],[],[],[],lb,ub,[],options);
[Rho_lam_opt,fval,exitflag,outpt] = Rho_lambda_opt(obs_x,theta,...
    obs_y,cntrl_input,calib_input,output,options)

%% Play with different stabilizing constant values
lambda = linspace(0.001,4,12);
for ii = 1:length(lambda)
    Sigma_eta = gp_cov(omega, D_in_x, D_in_x, rho, ...
        D_in_t, D_in_t, lambda(ii));

    % Get Sigma_D
    Sigma_D = Sigma_eta + padarray(Sigma_y,[m m],0,'post');
    % Add a nugget for computability
    WN = eye(size(Sigma_D)) * 10^(-4); Sigma_D = Sigma_D + WN;
    rcond(Sigma_D);

    % Get some needed values
    % DEPREC 
    c=10^3.6;
    lambda(ii)
    det_Sigma_D = det(c*Sigma_D)

    Sigma_D_inv = inv(Sigma_D);

    aa=D' * Sigma_D_inv * D
    det_S_D_12 = det_Sigma_D^(-1/2)
    expaa = exp(-1/2 * aa)
    Lik = det_S_D_12 * expaa
end
