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

%% Take a look at surfaces of example function output
theta1=linspace(0,3);
theta2=linspace(0,6);
% We need to convert these to a two-col matrix of all possible combinations
[ T1 T2 ] = meshgrid(theta1,theta2); 
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

%% Get simulation observations
n_cval  = 10 ; % Number of distinct c values to use
n_theta1 = 10 ; % Number of distinct theta1 values to use
n_theta2 = 10 ; % Number of distinct theta2 values to use
cvals  = linspace(1.5,2.5,n_cval)  ; % Get distinct c values
theta1vals = linspace(0,3,n_theta1) ; % Get distinct theta1 values
theta2vals = linspace(0,6,n_theta2) ; % Get distinct theta2 values

% Make a data frame sim_xt = [c theta1 theta2]
[ sc, st1, st2 ] = ndgrid(cvals,theta1vals,theta2vals) ;
sim_xt = [ sc(:) st1(:) st2(:) ] ;

% Get output
sim_dat = Ex_sim(sim_xt);

%% Run MCMC

% User defined values
M = 1e2;
desired_obs = [0 0 0];
which_outputs = [ 1 1 1 ] ; % Which of oscl, perf, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;
settings.Cost_lambda = 0; % remove prior on vf, thk
