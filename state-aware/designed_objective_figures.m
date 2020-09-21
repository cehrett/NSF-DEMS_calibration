% Designed objective figures
%% Set path string and add paths
clc; clear all; close all;

direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'state-aware']);

% Change dir
cd(dpath);


%% 3D view of univariate objective function
clc ; clearvars -except dpath ; close all ;

% Define minima and ranges of inputs
input_calib_min = 0 ; 
input_calib_range = 4 ;
input_cntrl_min = 2 ; 
input_cntrl_range = 1 ;
% Code for Univar_objective() automatically standardizes the output. We
% don't want that right now, so feed it fake output mean and std values.
output_mean = 0; 
output_sd   = 1;

% Get grids over inputs
M=2e1;
theta_grid = linspace(0,1,M);
x_grid = linspace(0,1,M);
[X,Theta] = meshgrid(x_grid,theta_grid);

% Get output
f=Univar_objective(X(:),Theta(:),0,...
    input_cntrl_min,input_cntrl_range,...
    input_calib_min,input_calib_range,...
    output_mean,output_sd);
f = reshape(f,M,M);

% Get surface
fig=figure();
srf = surf(X*input_cntrl_range + input_cntrl_min,...
    Theta*input_calib_range + input_calib_min,f,'EdgeAlpha',0.8);
ax=gca;
ax.View = [118.5000   15.0667];

% Save figure
% saveas(fig,'FIG_univ_SA_calib_prob_surf.png')
set(fig,'color','white');
% export_fig FIG_univ_SA_calib_prob_surf -png -m3 -painters 

%% Compare true optimum to mean result in univariate case
clc ; clearvars -except dpath ; close all ;

% Load the results
load([dpath,'state-aware\state-aware_results\'...
    '2018-11-29_sa_univariate_system'],...
    'results');

% Extract results
burn_in = results.settings.burn_in;
sampsrng = results.settings.input_calib_range;
sampsmin = results.settings.input_calib_min;
samps = results.theta1(burn_in+1:end,:) * sampsrng + sampsmin;
res = mean(samps);
reserr = std(samps);
truth=4/3 * (linspace(2,3,8) - 1);

% Make figure
fig = figure();
plot(linspace(2,3,8),truth,'r','LineWidth',2) ;
hold on ; 
errorbar(linspace(2,3,8),res,reserr,'o-','LineWidth',1.5); 
ylim([min(res-1.1*reserr) max(res+1.1*reserr)]);
xlabel('Control input (x)');
ylabel('Calibration input (\theta)');
lg = legend('True optima','Posterior mean, 1 s.d.','Location','northwest');

% Save figure
% saveas(fig,'FIG_univ_SA_calib_prob_surf.png')
set(fig,'color','white');
% export_fig FIG_univ_SA_results_vs_truth -png -m3 -painters 


%% Posterior marginal distributions of SA input at various controls
clc ; clearvars -except dpath ; close all ;

% Load the results
load([dpath,'state-aware\state-aware_results\'...
    '2018-11-29_sa_univariate_system'],...
    'results');

% Extract results
burn_in = results.settings.burn_in;
cal_rng = results.settings.input_calib_range;
cal_min = results.settings.input_calib_min;
samps = results.theta1(burn_in:end,:) * cal_rng + cal_min;
meantheta = mean(samps);
xlocs = linspace(2,3,8);
truetheta = 4/3 * (xlocs-1);

% Also get skew
skwnss = skewness(samps);


% Loop through each control setting, make a subplot each time
fig = figure('Pos',[10 10 800 450]);
for ii=1:size(samps,2)
    hax(ii) = subplot(2,4,ii);
    [f,xi]=ksdensity(samps(:,ii));
    plot(xi,f,'LineWidth',1.5);
    xlim([0.5,4]);ylim([0,1.75]);
    hold on;
    mp=plot([meantheta(ii) meantheta(ii)],...
        get(gca,'YLim'),'g:','LineWidth',1.5);
    tp=plot([truetheta(ii) truetheta(ii)],...
        get(gca,'YLim'),'r--','LineWidth',1);
    title(sprintf('x = %1.2f',xlocs(ii)));
    xlabel('\theta');
    txt = {'Skew:',sprintf('  %0.2f',skwnss(ii))};
    if ii==1 text(2.9,1.1,txt); else text(2.9,1.5,txt); end
end
suptitle(['Calibration parameter \theta posterior, '...
    'at eight control input locations x']);

%subplot(2,4,4);
lg = legend([mp,tp],'Mean','Truth','Location','southeast');
flushLegend(lg,hax(1),'northeast');


% Save figure
set(fig,'color','white');
% export_fig FIG_univ_SA_post_mean_marginals -png -painters -m3
