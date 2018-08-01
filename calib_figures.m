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
addpath([dpath,'Example\Ex_results']);

%% Show selection of desired observation near "elbow", and resulting calib.
clc ; clearvars -except dpath ; close all ;

%%% Load the preliminary CDO
load([dpath,'stored_data\'...
    '2018-07-25_discrepancy_d0'],...
    'results');
eouts = results.model_output.by_sample_est;

%%% Estimate the PF
[PF_os, PFidx] = nondominated(results.model_output.by_sample_est) ; 

%%% Put PF on standardized scale
omeans = mean(results.settings.output_means');
osds   = mean(results.settings.output_sds'  );
PF     = (PF_os - omeans)./ osds             ;

%%% Find closet point to des_obs
%orig_des_obs = results.settings.desired_obs  ;
orig_des_obs = [.74 .089 100 ] ; % This pt chosen to get at observed elbow
des_obs = (orig_des_obs - omeans)./osds      ;
[m,i] = min( sum( ( PF - des_obs ).^2, 2 ) ) ;
PF_optim = PF(i,:)                           ;

%%% Get new desired obs specified distance from PF in same dir as original
spec_dist = .2                               ;
dirvec_nonnormed = PF_optim - des_obs        ;
dirvec = dirvec_nonnormed/norm(dirvec_nonnormed) ;
des_obs_new = PF_optim - spec_dist * dirvec  ;
des_obs_new_os = des_obs_new .* osds + omeans; 

%%% Take a look
h=figure();
sc=scatter3(eouts(:,1),eouts(:,2),eouts(:,3),'g','MarkerEdgeAlpha',.2);
hold on;
scatter3(PF_os(:,1),PF_os(:,2),PF_os(:,3),'b','MarkerFaceColor','b',...
    'MarkerEdgeAlpha',.4,'MarkerFaceAlpha',.1)   ;                                ;
% scatter3(orig_des_obs(1),orig_des_obs(2),orig_des_obs(3))          ;
% line([orig_des_obs(1) PF_os(i,1)], [orig_des_obs(2) PF_os(i,2)], ...
%     [orig_des_obs(3) PF_os(i,3)])                                  ;
scatter3(des_obs_new_os(1),des_obs_new_os(2),des_obs_new_os(3),'r',...
    'MarkerFaceColor','r');
h.CurrentAxes.View = [-5.0000    5.2000];% [ 63 10] ;%[-8.4333 17.7333] ; 
title('Estimated Pareto front with desired observation');
xlabel('Deflection');ylabel('Rotation');zlabel('Cost');
set(h,'Color','w');
export_fig 'FIG_est_PF_with_des_obs' -png -m3
% saveas(h,'FIG_est_PF_with_des_obs.png');

%%% Now load the resulting calibration
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');
samps = results.samples_os(results.settings.burn_in+2:end,:) ;
outs =results.model_output.by_sample_est(results.settings.burn_in+2:end,:);
hh = copyobj(h,0); hold on;
scatter3(outs(:,1),outs(:,2),outs(:,3),30,'.m',...
    'MarkerFaceColor','m');
title('Posterior predictive distribution');

% saveas(hh,'FIG_post_pred_dist_with_model_range_and_des_obs.png');
export_fig 'FIG_post_pred_dist_with_model_range_and_des_obs' -png -m3

%% Get posterior scatterhist from calibration
clc ; clearvars -except dpath ; close all ;

%%% Load the calibration results
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');
samps = results.samples_os(results.settings.burn_in+2:end,:) ;

h=figure('pos',[10 10 410 310]);
sc=scatterhist(samps(:,1),samps(:,2),'Marker','.','MarkerSize',3);
xlim([0.2 0.6]); ylim([10 25]);
xlabel('Volume fraction');ylabel('Thickness');
title('Posterior distribution on \theta');
set(h,'Color','w');

% saveas(h,'FIG_post_dist_scatterhist.png');
export_fig 'FIG_post_dist_scatterhist' -png -m3;

%% Make post predictive plotmatrix
clc ; clearvars -except dpath ; close all ;

%%% Load the results
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');
outs = ...
    results.model_output.by_sample_est(results.settings.burn_in+2:end,:);
h=figure();
pm=plotmatrix(outs);
title('Predictive posterior distribution matrix');

%%% Throw on some labels
xsep=.05; ysep=.1;
ax1=h.Children(1);
xpos1 = ax1.XLim(1) + xsep * (ax1.XLim(2)-ax1.XLim(1)) ;
ypos1 = ax1.YLim(2) - ysep * (ax1.YLim(2)-ax1.YLim(1)) ;
ax2=h.Children(2);
xpos2 = ax2.XLim(1) + xsep * (ax2.XLim(2)-ax2.XLim(1)) ;
ypos2 = ax2.YLim(2) - ysep * (ax2.YLim(2)-ax2.YLim(1)) ;
ax3=h.Children(3);
xpos3 = ax3.XLim(1) + xsep * (ax3.XLim(2)-ax3.XLim(1)) ;
ypos3 = ax3.YLim(2) - ysep * (ax3.YLim(2)-ax3.YLim(1)) ;
text(ax1,xpos1,ypos1,'Deflection');
text(ax2,xpos2,ypos2,'Rotation');
text(ax3,xpos3,ypos3,'Cost');

%%% Save
% saveas(h,'FIG_pred_post_dist_plot_matrix.png');
set(h,'Color','w');
export_fig 'FIG_pred_post_dist_plot_matrix' -png -m3;

% text(.5,100,'test')

%% Get marginal plots from calibration with priors shown also
clc ; clearvars -except dpath ; close all ;

%%% Load the results
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');
samps = results.samples_os(results.settings.burn_in+2:end,:) ;

%%% Get the marginal plots
h1 = figure('rend','painters','pos',[10 10 310 210]) ; 
histogram(samps(:,1)) ;
% [f,xi]= ksdensity(samps(:,1),'Bandwidth',.005);
% plot(xi,f);