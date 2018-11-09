%% Code to reproduce figure 2 in paper.
% This code will produce a figure like that shown as figure 2 in the
% submitted paper. The figure is a scatterhist of the posterior
% distribution after calibration of the calibration parameters in the
% example case using simulated data. A code surrogate is not used here. As
% a scatterhist, the figure shows a scatter plot of the posterior
% distribution, along with marginal histograms along each axis. The figure
% also includes a contour plot of the distance, at each point in the
% parameter space, of the model output at that point from the target
% "desired observation" of [0 0 0].
%
% In order to produce this figure, the code here will complete calibration
% using the example model, with target desired observation 
% [0.7130 0.7144 17.9220]. 

addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\Technometrics\');
clc ; clearvars ; close all; 

% desired_obs is found by first choosing [0 0 0] as the performance target,
% and then updating along the line connecting [0 0 0] to the preliminary
% estimate of the pareto front so that the performance target lies closer
% to the model range. Here, the final performance target is specified.
desired_obs = [0.7130 0.7144 17.9220] ; 

% In the following line, we get settings for the MCMC routine. The inputs
% to MCMC_settings are the performance target, the number of control
% inputs, the number of calibration inputs, and the number of model inputs.
% Finally, we set the input minimima and ranges and the output means and 
% standard deviations, so that the inputs and outputs can be normalized and
% standardized (respectively).
settings = MCMC_settings(desired_obs,1,2,3,...
    'input_cntrl_mins',1.9500,...
    'input_calib_mins',[0 0],...
    'input_cntrl_ranges',0.1,...
    'input_calib_ranges',[3 6],...
    'output_means',[0.9361 ; 0.8107 ; 21.1364],...
    'output_sds',[0.0605 ; 0.1053 ; 3.4961]);

% Now we perform the calibration.
results = MCMC_calib(settings) ; 

%% Now make the figure
% The figure sets the posterior distribution against a contour plot showing
% the distance, at each point in the parameter space [0,3]\times[0,6], of
% the average model output at that point from the desired model output.
% Thus, we hope that our posterior distribution peaks at the point in the
% parameter space that is closest to the desired output, and that more
% generally the posterior distribution describes the region of the
% parameter space that is "nearest" to the desired observation.
%
% To make the figure, then, we must determine the distance of the average
% model output from the desired observation at each point of the parameter
% space. To do this, since this is a simulation involving a relatively
% computationally cheap model, we simply run the model at each point of a
% large grid over the parameter space, and use the resulting output to draw
% the contour plot. The results of that grid evaluation of the model have
% already been completed and are here loaded from stored data.
% Load true samples;
load('direct_fn_output');

% Put the outputs and the desired observation on the standardized scale
meanout = results.settings.output_means' ;
sdout   = results.settings.output_sds' ;
cost_std = (ctheta_output(:,6) - meanout(3))/...
    sdout(3);
defl_std = (ctheta_output(:,4) - meanout(1))/...
    sdout(1);
rotn_std = (ctheta_output(:,5) - meanout(2))/...
    sdout(2);
outputs_std = [defl_std rotn_std cost_std];
desired_obs = (desired_obs-meanout) ./ sdout;

% True optimum: find using direct function output
zero_stdized = ([0 0 0] - meanout) ./ sdout;
zero_dists = sqrt( sum( (outputs_std - zero_stdized).^2 , 2) );
[m,idx] = min(zero_dists);
optim = ctheta_output(idx,2:3);

% Gather posterior samples
samps = results.samples_os(results.settings.burn_in+2:end,:);

f=figure();

% Now get Euclidean norms of each standardized output
dists = sqrt ( sum ( (outputs_std-desired_obs).^2 , 2 ) ) ;
redists = reshape(dists,1000,1000);
theta1 = reshape(ctheta_output(:,2),1000,1000);
theta2 = reshape(ctheta_output(:,3),1000,1000);

%%% Make scatterhist of posterior samples
colormap autumn
sc=scatterhist(samps(:,1),samps(:,2),'Marker','.','Color','b',...
    'Markersize',1); 
hold on;
title(['Posterior \theta samples: desired obs. '...
    '[0.71 0.71 17.92], \lambda_\delta = 1']);

% Now add contour plot 
%[theta1,theta2] = meshgrid(linspace(0,3,1000),linspace(0,6,1000));
[C,h]= contour(theta1,theta2,redists,[1 2 3 4 5 ],'LineWidth',3);
clabel(C,h,'fontsize',12);
xlabel('\theta_1'); ylabel('\theta_2');

%%% Add true optimum
p=plot(optim(1),optim(2),'ok','MarkerSize',7,'MarkerFaceColor','m',...
    'LineWidth',2);


