% Optimize omega, rho wrt average lambda

%% Get data
load raw_dat.mat
tdat = Tdat(raw_dat,3); % rescaling inputs, standardizing outputs

%% Covariance parameter settings
omega = [0.2 0.82];
rho   = [0.39 0.99];
lambda = 0.5;

%% Assign variables
theta = rho;
obs_x = [];
obs_y = [];
sim_x = tdat.input(:,1:2);
sim_t = tdat.input(:,3:4);
sim_y = tdat.output;
cntrl_input = sim_x;
calib_input = sim_t;
output = sim_y;
c     = 10^3.6;

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
    calib_input,output,c),1.01,1.015);

%% Define optimization fn of Rho over the integral (given data, theta, etc)
Rho_opt = @(obs_x,theta,obs_y,cntrl_input,calib_input,output,c) ...
    fminsearch( @(Rho) -intgrl(Rho,obs_x,theta,obs_y,cntrl_input, ...
    calib_input,output,c),[omega' ; rho' ]);

%% Perform optimization
Rho_optimum = Rho_opt(obs_x,theta,obs_y,cntrl_input,...
    calib_input,output,c)

save('Rho_optimum','Rho_optimum');
