%% Set path string and add paths
clc; clear all; close all;

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'state-aware']);

%% MCMC chains of state-aware variable t1
clc ; clearvars -except dpath ; close all ;

% Load results
load([dpath,'state-aware\state-aware_results\'...
    '2018-12-08_sa_true_fn_ld1_t1calib_x23'],...
    'results');

f = figure('pos',[10 10 600 480 ]);
hold on;

% Add colors
colors = [      0         0    1.0000 ;
    1.0000         0         0 ;
         0    1.0000         0 ;
         0         0    0.1724 ;
    1.0000    0.1034    0.7241 ;
    1.0000    0.8276         0 ;
         0    0.3448         0 ;
    0.5172    0.5172    1.0000];
% [ 1    0    0 ; 
%     0.5767    0.7127    0.0967 ;
%     0.1829    0.5005    0.8181 ;
%     0.5399    0.2711    0.8175 ;
%     0.8865    0.0596    0.7224 ;
%     0.3000    1         0.1500 ;
%     0.0899    0.6424    0.9596 ;
%     0.0800    0.0800    0.3000 ];

% Plot each chain
for ii = 1:8
    plot( results.theta1(results.settings.burn_in+1:end,ii),...
        'LineWidth',1.25,'Color',colors(ii,:));
    lgcell{ii} = num2str(...
        results.settings.cntrl_input(ii) * ...
        results.settings.input_cntrl_min + ...
        results.settings.input_cntrl_range, 'x=%3.2f');
end
lg1 = legend(lgcell);
title(...
    'MCMC samples of state-aware \theta_1 at each of 8 control inputs x');

% Save figure
set(f,'Color','w');
% export_fig FIG_state-aware_t1 -png -m3 -painters 

% Now zoom in
% xlim([11000 12000]);

% Save figure
set(f,'Color','w');
% export_fig FIG_state-aware_t1_zoom -png -m3 -painters 

%% MCMC chains of state-aware variable t2
clc ; clearvars -except dpath ; close all ;

% Load results
load([dpath,'state-aware\state-aware_results\'...
    '2018-10-19_sa_true_fn_ld1_t2calib'],...
    'results');

f = figure('pos',[10 10 600 480 ]);
hold on;

% Add colors
colors = [ 1    0    0 ; 
    0.5767    0.7127    0.0967 ;
    0.1829    0.5005    0.8181 ;
    0.5399    0.2711    0.8175 ;
    0.8865    0.0596    0.7224 ;
    0.3000    1         0.1500 ;
    0.0899    0.6424    0.9596 ;
    0.0800    0.0800    0.3000 ];

% Plot each chain
for ii = 1:8
    plot( results.theta1(results.settings.burn_in+1:end,ii),...
        'LineWidth',1.25,'Color',colors(ii,:));
    lgcell{ii} = num2str(...
        results.settings.cntrl_input(ii) * ...
        results.settings.input_cntrl_min + ...
        results.settings.input_cntrl_range, 'x=%3.2f');
end
lg1 = legend(lgcell);
title(...
    'MCMC samples of state-aware \theta_2 at each of 8 control inputs x');

% Save figure
set(f,'Color','w');
% export_fig FIG_state-aware_t2 -png -m3 -painters 

% Now zoom in
xlim([9000 10000]);

% Save figure
set(f,'Color','w');
export_fig FIG_state-aware_t2_zoom -png -m3 -painters 
