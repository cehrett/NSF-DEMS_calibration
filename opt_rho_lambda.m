function orl = opt_rho_lambda(sim_xt,sim_dat,num_cal,omega,rho,lambda)
% This function takes simulation input and corresponding output, and
% returns hyperparameters rho and lambda for a GP emulator for that data,
% following the framework of Willams et al 2006
% sim_xt is expected to be in the form of control variable columns followed
% by calibration variable columns. num_cal gives the number of
% calibration variables.
% omega, rho and lambda are initial values for the optimization.


%% Normalize the inputs and standardize the outputs
tdat = Tdat([sim_xt sim_dat],size(sim_dat,2));
% Divide up tdat. Note that the code currently does not support inclusion
% of field observations, so we set those to be empty.
obs_x = [];
obs_y = [];
theta = .5 * ones(1,num_cal); % This is a placeholder, the value doesn't 
                              % matter because, again, the code currently
                              % doesn't support field observations
num_cntrl = size(tdat.input,2) - num_cal ; % Infer # of cntrl vars
sim_x = tdat.input(:,1:num_cntrl);
sim_t = tdat.input(:,num_cntrl+1:num_cntrl+num_cal);
sim_y = tdat.output;

%% Perform the optimization
% Set fmincon options
options = optimset('Display','iter','PlotFcns',@optimplotfval,...
    'TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',2000);
% Set initial value
Rho_lam_init = [omega rho lambda]; 
% Set upper and lower boundaries. Each rho has [0,1] support; lambda has
% positive support.
lb = zeros(1,num_cntrl+num_cal+1);
ub = lb + 1 ; ub(length(ub)) = Inf;

% Define the optimizing function (in terms of fmincon)
Rho_lambda_opt = @(obs_x,theta,obs_y,cntrl_input,...
    calib_input,output,options)...
    fmincon(@(Rho_lam) -logLorl(Rho_lam(1:num_cntrl),...
    Rho_lam(num_cntrl+1:num_cntrl+num_cal),...
    Rho_lam(num_cntrl+num_cal+1),...
    obs_x,theta,obs_y,cntrl_input,calib_input,output),...
    Rho_lam_init,[],[],[],[],lb,ub,[],options);

% Call the optimizing function
[Rho_lam_opt,~,~,~] = Rho_lambda_opt(obs_x,theta,...
    obs_y,sim_x,sim_t,sim_y,options)

orl = Rho_lam_opt;

end