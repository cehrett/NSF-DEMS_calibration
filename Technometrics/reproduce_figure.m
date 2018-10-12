% Code to reproduce figure 2 in paper.
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
% using the example model, with target desired observation [0 0 0]. 

clc ; clearvars -except dpath ; close all; 

% Get settings
desired_obs = [0 0 0 ] ; %[ 0.7130 0.7144 17.9220 ] ; 
settings = MCMC_settings(desired_obs,[nan],[nan nan],[nan nan nan],...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant');





