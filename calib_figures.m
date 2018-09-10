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
n=700;
burn_in = results.settings.burn_in;
eouts = results.model_output.by_sample_est(burn_in:burn_in+n,:);

%%% Estimate the PF
[PF_os, PFidx] = nondominated(eouts) ; 
PFidx = PFidx + results.settings.burn_in; % get idx in original full set

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
sc=scatter3(eouts(:,1),eouts(:,2),eouts(:,3),'g','MarkerEdgeAlpha',1,...
    'MarkerFaceAlpha',.2,'MarkerFaceColor','g');
hold on;
scatter3(PF_os(:,1),PF_os(:,2),PF_os(:,3),'b','MarkerFaceColor','b',...
    'MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.2)   ;
% scatter3(orig_des_obs(1),orig_des_obs(2),orig_des_obs(3))          ;
% line([orig_des_obs(1) PF_os(i,1)], [orig_des_obs(2) PF_os(i,2)], ...
%     [orig_des_obs(3) PF_os(i,3)])                                  ;
scatter3(des_obs_new_os(1),des_obs_new_os(2),des_obs_new_os(3),'r',...
    'MarkerFaceColor','r');
h.CurrentAxes.View = [-3.9333   10.5333] ; 
% [-5.0000    5.2000];% [ 63 10] ;%[-8.4333 17.7333] ; 
title('Estimated Pareto front with desired observation');
xlabel('Deflection');ylabel('Rotation');zlabel('Cost');
set(h,'Color','w');
% export_fig 'FIG_est_PF_with_des_obs' -png -m3
% saveas(h,'FIG_est_PF_with_des_obs.png');

%%% Now load the resulting calibration
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');
samps = results.samples_os(results.settings.burn_in+2:end,:) ;
outs =results.model_output.by_sample_est(results.settings.burn_in+2:end,:);
hh = copyobj(h,0); hold on;
sc1=scatter3(outs(:,1),outs(:,2),outs(:,3),30,'.m',...
    'MarkerFaceColor','m');
title('Posterior predictive distribution');

% saveas(hh,'FIG_post_pred_dist_with_model_range_and_des_obs.png');
% export_fig 'FIG_post_pred_dist_with_model_range_and_des_obs' -png -m3

%%% Now color the closest alpha percent of the posterior dist to des obs
alpha=0.5;
souts = (outs - omeans)./osds;
posdiffs = max(0, souts-des_obs_new);
dists = sum( (posdiffs).^2, 2);
close_idx = dists < quantile(dists,alpha);
hhh = copyobj(h,0); hold on;
scatter3(outs(close_idx,1),outs(close_idx,2),outs(close_idx,3),30,'.b',...
    'MarkerFaceColor','b');
scatter3(outs(~close_idx,1),outs(~close_idx,2),outs(~close_idx,3),30,'.m',...
    'MarkerFaceColor','m');

% %%% Add code to turn the first figure in this section into a version for
% %%% the SCSC poster
% h.Position = [360.3333  197.6667  452.0000  259.3333];
% set(h,'Color',[251/255 244/255 245/255]);

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
h1 = figure('rend','painters','pos',[10 10 610 160]) ; 
subplot(1,2,1);
histogram(samps(:,1), 'Normalization','pdf') ;
xlim([0.2 0.6]);
unifval = 1/.4;
hold on;
plot([0.2 0.6], [unifval unifval],'--r','LineWidth',2);
title('Volume fraction');

subplot(1,2,2);
histogram(samps(:,2), 'Normalization','pdf') ;
xlim([10 25]);
unifval = 1/15;
hold on;
plot([10 25], [unifval unifval],'--r','LineWidth',2);
title('Thickness (mm)');

%%% Save
set(h1,'Color','none');
export_fig FIG_posterior_marginals_with_priors -png -m3;

%% Get Pareto bands figure from cost_grid analysis using discrepancy
clc ; clearvars -except dpath ; close all ; 

%%% Load the cost_grid results
load([dpath,'stored_data\'...
    '2018-08-03_cost_grid_discrepancy_results'],...
    'results');

% Collect Cost_lambdas, and posterior mean and sds for costs, defl, rot, as
% well as upper and lower .05 quantiles
m=size(results,1); % Store number of target cost_lambdas
cost_lambda = zeros(m,1); % This will store cost_lambdas
cred_level = 90; % Set desired level for credible bands (in %)
alpha = (100-cred_level)/100; % Convert cred_level to alpha level
pmo = zeros(m,3); % This will store posterior mean output of emulator
pdo = zeros(m,3); % ``'' median output
pso = zeros(m,3); % ``'' appropriate multiple of standard deviations
plo = zeros(m,3); % ``'' lower (alpha/2) quantile
puo = zeros(m,3); % ``'' upper (alpha/2) quantile
for ii = 1:m % This loop populates the above arrays
    pmo(ii,:) = results{ii}.post_mean_out;
    pdo(ii,:) = quantile(results{ii}.model_output.by_sample_est,0.5);
    pso(ii,:) = norminv(1-alpha/2) * ...
        mean(results{ii}.model_output.by_sample_sds);
    plo(ii,:) = quantile(results{ii}.model_output.by_sample_est,alpha/2);
    puo(ii,:) = quantile(results{ii}.model_output.by_sample_est,1-alpha/2);
    cost(ii) = results{ii}.desired_obs(3);
end
% Now we break the arrays up each into 3 vectors, one for each output
post_cost_mean = pmo(:,3);
post_defl_mean = pmo(:,1);
post_rotn_mean = pmo(:,2);
post_cost_median = pdo(:,3);
post_defl_median = pdo(:,1);
post_rotn_median = pdo(:,2);
post_cost_sd = pso(:,3);
post_defl_sd = pso(:,1);
post_rotn_sd = pso(:,2);
post_cost_lq = plo(:,3);
post_cost_uq = puo(:,3);
post_defl_lq = plo(:,1);
post_defl_uq = puo(:,1);
post_rotn_lq = plo(:,2);
post_rotn_uq = puo(:,2);
% Get quantiles plus code uncertainty
post_cost_uq_cu = post_cost_uq + post_cost_sd;
post_cost_lq_cu = post_cost_lq - post_cost_sd;
post_defl_uq_cu = post_defl_uq + post_defl_sd;
post_defl_lq_cu = post_defl_lq - post_defl_sd;
post_rotn_uq_cu = post_rotn_uq + post_rotn_sd;
post_rotn_lq_cu = post_rotn_lq - post_rotn_sd;
% Get ylims for the two sets of plots
ylimrat=1.01;
ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];
ylim_rotn = [min(post_rotn_lq_cu)/ylimrat max(post_rotn_uq_cu)*ylimrat];

% Begin figures
h=figure('rend','painters','pos',[10 10 800 400]);
x = 96:1:350; % x fills the cost domain
% Now begin plot 1/3
subplot(1,2,2)
% Get main curve
pdefl = pchip(cost,post_defl_mean,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(cost,post_defl_uq,x);
pdefllq = pchip(cost,post_defl_lq,x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(f,'facealpha',.25,'EdgeAlpha',.25);
hold on;
plot(...%cost,post_defl_mean,'or',...
    x,pdefl,'-r',...
    x,pdefluq,':k',...
    x,pdefllq,':k');
xl2=xlabel('Target cost');
ylabel('Deflection');
xlim([96,350]);
ylim(ylim_defl);

% % Here's plot 2/3
% subplot(1,3,2)
% % Get main curve
% protn = pchip(cost,post_rotn_mean,x);
% % Get upper and lower 0.05 quantiles curves
% protnuq = pchip(cost,post_rotn_uq,x);
% protnlq = pchip(cost,post_rotn_lq,x);
% f=fill([ x , fliplr(x) ], [protnuq, fliplr(protnlq)],'k');
% set(f,'facealpha',.25);
% hold on;
% plot(cost,post_rotn_mean,'or',...
%     x,protn,'-r',...
%     x,protnuq,'-k',...
%     x,protnlq,'-k');
% xl3=xlabel('Target cost');
% ylabel('Rotation');
% xlim([96,350]);
% ylim(ylim_rotn);

% Here's plot 3/3
subplot(1,2,1)
% Get main curve
pcost = pchip(cost,post_cost_mean,x);
% Get upper and lower 0.05 quantiles curves
pcostuq = pchip(cost,post_cost_uq,x);
pcostlq = pchip(cost,post_cost_lq,x);
go_fill_unc=fill([ x , fliplr(x) ], [pcostuq, fliplr(pcostlq)],'k');
set(go_fill_unc,'facealpha',.25,'edgealpha',.25);
hold on;
go_plot_mean=plot(x,pcost,'-r');
plot(...%cost,post_cost_mean,'or',...%x,pcost,'-r',...
    x,pcostuq,':k',...
    x,pcostlq,':k');
% Plot 2sd errbar
% errorbar(cost_lambda,post_cost_mean,post_cost_sd,'ob'); 
xl1=xlabel('Target cost');
ylabel('Observed cost');
xlim([96,350]);
ylim(ylim_cost);
%plot(x,x,'-k','LineWidth',2);

% % Save the figure temporarily so we can mess with it later, because
% % suptitle seems to mess things up somehow for making changes after calling
% % it
% savefig(h,'tempfig');
% 
% % Now add a main title and fix any infelicities
% suptitle(['Deflection vs. (known) target cost,',...
%     ' with ',num2str(cred_level),'% credible interval']); 
% p = get(xl1,'position');
% set(xl1,'position',p + [0 2.75 0]);
% p = get(xl2,'position');
% set(xl2,'position',p + [0 0.00125 0])
% % p = get(xl3,'position');
% % set(xl3,'position',p + [0 0.0002 0])
figpos = get(h,'pos');

% saveas(h,'FIG_cost_grid_pareto.png');

% Now add in code uncertainty. That is, the above assumes that the GP
% emulator nails the FE code precisely. But of course the GP emulator has
% nonnegligible variance. That's the code uncertainty. So our confidence
% bands should reflect it. So we add it in here, by dropping the
% appropriate multiple of the sd from each lower quantile and adding it to
% each upper quantile.
% First, open the figure prior to calling suptitle.
% h=openfig('tempfig');
subplot(1,2,2);
pdefluq_code_uncert = pchip(cost,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(cost,post_defl_lq_cu,x);
f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
    fliplr(pdefluq)],'b');
ff=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
    fliplr(pdefllq)],'b');
set(f,'facealpha',.25,'EdgeAlpha',.25);
set(ff,'facealpha',.25,'EdgeAlpha',.25);
xl2=xlabel('Target cost');
ylim(ylim_defl);

% subplot(1,3,2);
% protnuq_code_uncert = pchip(cost ,post_rotn_uq_cu,x);
% protnlq_code_uncert = pchip(cost ,post_rotn_lq_cu,x);
% f=fill([ x , fliplr(x) ], [protnuq_code_uncert,...
%     fliplr(protnuq)],'b');
% ff=fill([ x , fliplr(x) ], [protnlq_code_uncert,...
%     fliplr(protnlq)],'b');
% set(f,'facealpha',.25);
% set(ff,'facealpha',.25);
% xl3=xlabel('Target cost');
% ylim(ylim_rotn);

subplot(1,2,1);
pcostuq_code_uncert = pchip(cost ,post_cost_uq_cu,x);
pcostlq_code_uncert = pchip(cost ,post_cost_lq_cu,x);
go_fill_cunc_up=fill([ x , fliplr(x) ], [pcostuq_code_uncert,...
    fliplr(pcostuq)],'b');
go_fill_cunc_dn=fill([ x , fliplr(x) ], [pcostlq_code_uncert,...
    fliplr(pcostlq)],'b');
set(go_fill_cunc_up,'facealpha',.25,'EdgeAlpha',.25);
set(go_fill_cunc_dn,'facealpha',.25,'EdgeAlpha',.25);
xl1=xlabel('Target cost');
ylim(ylim_cost);
go_plot_diag=plot(ylim_cost,ylim_cost,'--b','LineWidth',2);

% Now add a main title and fix any infelicities
suptitle(['Posterior estimate vs. target cost,',...
    ' with ',num2str(cred_level),'% credible interval ']); 
set(h,'pos',figpos); % Just so we can reuse the positioning code from above
p = get(xl1,'position');
set(xl1,'position',p + [0 2.75 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.00125 0])
% p = get(xl3,'position');
% set(xl3,'position',p + [0 0.0002 0])

% Now add a legend.
leg_gos = [go_plot_mean go_fill_unc go_fill_cunc_up go_plot_diag];
%leg_str = ['Posterior predictive mean'...
%     '90\% credible interval w/o code uncertainty'...
%     'Extension of 90\% credible interval for code uncertainty'...
%     'Diagonal for reference'];
lg=legend(leg_gos,'Posterior predictive mean',...
    'C.I. w/o code uncertainty',...
    sprintf('Expanded C.I. w/ code uncertainty'),...
    'Diagonal for reference',...
    'Location','northwest');
lg.Position(1:2)=[.612 .6875];

%%% Change to version for SCSC poster
h.Position = [ 20 20 510 320];
lg.String = {['Post. pred. mean'],...
    ['C.I. w/o code' newline 'uncertainty'],...
    ['Expanded C.I. w/' newline 'code uncertainty'],...
    'Diagonal for ref.'};
lg.Position = [0.635    0.597    0.2595    0.2089];
set(h,'Color','none');
export_fig FIG_cost_grid_pareto_bands -png -m3 -painters

    
%%% Save
% set(h,'Color','white');
%export_fig FIG_cost_grid_pareto_bands -png -m3 -painters
%saveas(h,'FIG_cost_grid_pareto_with_code_uncert.png');

%% Compare prior predictive distribution to posterior predictive distributn
clc ; clearvars -except dpath ; close all ;

%%% Load prior predictive results
load([dpath,'stored_data\'...
    '2018-09-03_prior_pred_distrib'],...
    'prior_pred_dist');
prsamps = prior_pred_dist.prior_pred_pts;
clear prior_pred_dist;

%%% Load calib results
load([dpath,'stored_data\'...
    '2018-07-27_discrepancy_d-elbow_d-p2'],...
    'results');
posamps = results.model_output.by_sample_est;
des_obs = results.settings.desired_obs;
clear results; 

%%% Make figure using histograms
f=figure('pos',[10 10  580.0000  246.6667]);
% Deflection
subplot(1,3,1);
histogram(posamps(:,1),'Normalization','pdf','Edgecolor','none'); 
hold on;
histogram(prsamps(:,1),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Deflection','VerticalAlignment','bottom','Units','normalized');
text(1.715,102,'Rotation','VerticalAlignment','bottom');
xlim([0.6 0.85]);
ylim([0 110]);
line([des_obs(1) des_obs(1)],ylim,'Color','black','Linestyle','--',...
    'linewidth',2);
% Rotation
subplot(1,3,2);
histogram(posamps(:,2),'Normalization','pdf','Edgecolor','none'); 
hold on;
histogram(prsamps(:,2),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Rotation','VerticalAlignment','bottom','Units','normalized');
xlim([0.075,0.105])
ylim([0 700]);
line([des_obs(2) des_obs(2)],ylim,'Color','black','Linestyle','--',...
    'linewidth',2);
% Cost
subplot(1,3,3);
histogram(posamps(:,3),'Normalization','pdf','Edgecolor','none'); 
hold on;
histogram(prsamps(:,3),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Cost','VerticalAlignment','bottom','Units','normalized');
ylim([0 .0700]);
line([des_obs(3) des_obs(3)],ylim,'Color','black','Linestyle','--',...
    'linewidth',2);
% Add suptitle
st=suptitle('Prior (red) and posterior (blue) predictive distributions');
st.Position=[0.5 -.1 0];

%%% Save
set(f, 'Color','none');
export_fig FIG_prior_vs_posterior_dist -png -m3 -painters