% Simulation workspace

clc; clear all; close all;

%% Set path string
direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

%% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'Example\Ex_results']);

%% Take a look at surfaces of example function output
theta1=linspace(0,3);
theta2=linspace(0,6);
% We need to convert these to a two-col matrix of all possible combinations
[ T1 , T2 ] = meshgrid(theta1,theta2); 
ths = cat(2,T1',T2');
Theta = reshape(ths,[],2);

% Set other parameters
c = repmat(1.5,size(Theta,1),1);

% Now get output values
opc = Ex_sim([c Theta]);

% Normalize the outputs
oscl_n = (opc(:,1)-min(opc(:,1))) / range(opc(:,1)) ;
perf_n = (opc(:,2)-min(opc(:,2))) / range(opc(:,2)) ;
cost_n = (opc(:,3)-min(opc(:,3))) / range(opc(:,3)) ;

% Now take a look at the surfaces
oscls = reshape(oscl_n,[],length(theta1));
perfs = reshape(perf_n,[],length(theta1));
costs = reshape(cost_n,[],length(theta1));

ec = 'black' ;  % edge color
ea = .25       ;  % edge alpha
fa = .75       ;  % face alpha

surf(theta2,theta1,oscls,'FaceColor','red','EdgeColor',ec,...
    'EdgeAlpha',ea,'FaceAlpha',fa);
axis vis3d;

hold on;
surf(theta2,theta1,perfs,'FaceColor','blue','EdgeColor',ec,...
    'EdgeAlpha',ea,'FaceAlpha',fa); 
axis vis3d;

surf(theta2,theta1,costs,'FaceColor','green','EdgeColor',ec,...
    'EdgeAlpha',ea,'FaceAlpha',fa); 
axis vis3d;
xlabel('\theta_2'); ylabel('\theta_1'); zlabel('Outcomes');

% Add line at posterior mean after below calibration
pmo = mean(samples(settings.burn_in:end,:)) .* settings.input_calib_ranges;
hold on;
plot3([pmo(1) pmo(1)], [pmo(2) pmo(2)], get(gca,'Zlim'), 'k',...
    'LineWidth',6);

%% Get simulation observations
n_cval  = 5 ; % Number of distinct c values to use
n_theta1 = 7 ; % Number of distinct theta1 values to use
n_theta2 = 7 ; % Number of distinct theta2 values to use
cvals  = linspace(1.5,2.5,n_cval)  ; % Get distinct c values
theta1vals = linspace(0,3,n_theta1) ; % Get distinct theta1 values
theta2vals = linspace(0,6,n_theta2) ; % Get distinct theta2 values

% Make a data frame sim_xt = [c theta1 theta2]
[ sc, st1, st2 ] = ndgrid(cvals,theta1vals,theta2vals) ;
sim_xt = [ sc(:) st1(:) st2(:) ] ;

% Get output
sim_y = Ex_sim(sim_xt);

%% Run MCMC

% User defined values
M = 1e3;
desired_obs = [0 0 0];
which_outputs = [ 1 1 1 ] ; % Which of oscl, perf, cost
% Calculated optimum for example simulation:
Rho_lam_optimum  = [  0.280981573480363   0.999189406633873...
   0.600440750045477  0.719652153362981   0.102809702497319...
   0.000837772517865 ] ;

% Settings
settings = MCMC_settings (M,desired_obs,sim_xt(:,1),sim_xt(:,2:3),...
    sim_y,which_outputs,Rho_lam_optimum);
settings.doplot = false;
settings.doplot = true;

% MCMC
[samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);

% Get extra info about results and save everything
post_mean_out = em_out(samples,settings)
results = struct('samples',samples,...
    'sigma2',sigma2_rec,...
    'Sigma',Sigma,...
    'init',samples(1,:),...
    'desired_obs',desired_obs,...
    'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
    'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
    'post_mean_out',post_mean_out,...
    'settings',settings);

save([dpath,'Example\Ex_results\'...
'2018-05-11_d0_incl_min_cost'],...
'results');


%% Gather results over grid of cost values
cost_grid = linspace(15,30,8);

for ii = 1:length(cost_grid)
    
        
    
end
