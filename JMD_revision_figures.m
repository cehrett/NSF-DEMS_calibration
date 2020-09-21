% Figures for the revision for the JMD submission of CTO paper.

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
addpath([dpath,'Example\Ex_results']);

% Change dir
cd(dpath);


%% WTA Pareto bands (new 2020-04-23)
clc ; clearvars -except dpath ; close all ;

%%% Load the results
locstr= ['C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\'...
    'stored_data\2020-04-26_CTO_costgrid_size30'];
load(locstr);
mean_y = results{1}.settings.mean_y ; std_y = results{1}.settings.std_y ;

% Collect Cost_lambdas, and posterior mean and sds for costs, defl, rot, as
% well as upper and lower .05 quantiles
m=length(results); % Store number of target cost_lambdas
cost_lambda = zeros(m,1); % This will store cost_lambdas
cred_level = 90; % Set desired level for credible bands (in %)
alpha = (100-cred_level)/100; % Convert cred_level to alpha level
pmo = zeros(m,2); % This will store posterior mean output of emulator
pdo = zeros(m,2); % ``'' median output
pso = zeros(m,2); % ``'' appropriate multiple of standard deviations
plo = zeros(m,2); % ``'' lower (alpha/2) quantile
puo = zeros(m,2); % ``'' upper (alpha/2) quantile
for ii = 1:m % This loop populates the above arrays
    output_means = results{ii}.model_output.means .* std_y + mean_y;
    output_sds = sqrt(results{ii}.model_output.vars) .* std_y;
    pmo(ii,:) = mean(output_means);
    pdo(ii,:) = quantile(output_means,0.5);
    pso(ii,:) = norminv(1-alpha/2) * ...
        mean(output_sds);
    plo(ii,:) = quantile(output_means,alpha/2);
    puo(ii,:) = quantile(output_means,1-alpha/2);
    cost(ii) = results{ii}.settings.obs_y(1,end) * std_y(2) + mean_y(2);
end
% Now we break the arrays up each into 2 vectors, one for each output
post_cost_mean = pmo(:,2);
post_defl_mean = pmo(:,1);
post_cost_median = pdo(:,2);
post_defl_median = pdo(:,1);
post_cost_sd = pso(:,2);
post_defl_sd = pso(:,1);
post_cost_lq = plo(:,2);
post_cost_uq = puo(:,2);
post_defl_lq = plo(:,1);
post_defl_uq = puo(:,1);
% Get quantiles plus code uncertainty
post_cost_uq_cu = post_cost_uq + post_cost_sd;
post_cost_lq_cu = post_cost_lq - post_cost_sd;
post_defl_uq_cu = post_defl_uq + post_defl_sd;
post_defl_lq_cu = post_defl_lq - post_defl_sd;
% Get ylims for the two sets of plots
ylimrat=1.01;
ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];

%%% Begin figures
% Set alphas for two types of uncertainty
alpha_wcu = 0.5;  %with code uncertainty
alpha_wocu= 0.15; %without
h=figure('rend','painters','pos',[10 10 360 240]);
x = 96:1:350; % x fills the cost domain
[subplts , pos] = tight_subplot(1,2,0.175,[0.15 0.02],[0.11 0.01]);

% Now begin plot 1/2
axes(subplts(1));
% Get main curve
pdefl = pchip(cost,post_defl_median,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(cost,post_defl_uq,x);
pdefllq = pchip(cost,post_defl_lq,x);
unc_wo_cu = fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(unc_wo_cu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
hold on;
median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
% plot(x,pdefluq,':k',...
%      x,pdefllq,':k');
xl2=xlabel('Target cost');
ylabel('Deflection');
xlim([96,350]);
ylim(ylim_defl);

figpos = get(h,'pos');

% Add NSGA-II results
locstr = [dpath,'stored_data\'...
    '2020-04-23_NSGA2_results'];
load(locstr);
nsga2_res = plot(result.final_obj(:,2),result.final_obj(:,1),'*');


% Now add a legend.
ylim(ylim+[0 0.03]);
leg_gos = [nsga2_res median_line ];% go_plot_diag];
lg=legend(leg_gos,sprintf('NSGA-II results'),...
    sprintf('Posterior\npredictive median'),...
    'Location','northeast');
% lg.Position(1:2)=[.623 .725];
flushLegend(lg,'northeast');
lg.Box='off';
% lgpos = lg.Position;
% lg.Position = lgpos + [-.004 -.002 -.004 -.002];


% saveas(h,'FIG_cost_grid_pareto.png');

% Now add in code uncertainty. That is, the above assumes that the GP
% emulator nails the FE code precisely. But of course the GP emulator has
% nonnegligible variance. That's the code uncertainty. So our confidence
% bands should reflect it. So we add it in here, by dropping the
% appropriate multiple of the sd from each lower quantile and adding it to
% each upper quantile.
% First, open the figure prior to calling suptitle.
% h=openfig('tempfig');
axes(subplts(2));
pdefluq_code_uncert = pchip(cost,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(cost,post_defl_lq_cu,x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(f,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
hold on;
f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
    fliplr(pdefluq)],'k');
unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
    fliplr(pdefllq)],'k');
hold on;
median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
set(f,'facealpha',alpha_wcu,'EdgeAlpha',alpha_wcu);
set(unc_w_cu,'facealpha',alpha_wcu,'EdgeAlpha',alpha_wcu);
ylabel('Deflection');
xl2=xlabel('Target cost');
xlim([183.7013 186.1534]);
ylim([0.7255 0.7273]);

% Add NSGA-II results
locstr = [dpath,'stored_data\'...
    '2020-04-23_NSGA2_results'];
load(locstr);
plot(result.final_obj(:,2),result.final_obj(:,1),'*');

% Now add a main title and fix any infelicities
% suptitle(['Posterior estimate vs. target cost,',...
%     ' with ',num2str(cred_level),'% credible interval ']); 
set(h,'pos',figpos); % Just so we can reuse the positioning code from above
% p = get(xl1,'position');
% set(xl1,'position',p + [0 2.75 0]);
% p = get(xl2,'position');
% set(xl2,'position',p + [0 0.00125 0])
% p = get(xl3,'position');
% set(xl3,'position',p + [0 0.0002 0])

% Now add a legend.
ylim(ylim+[0 0.0007]);
leg_gos = [median_line unc_wo_cu unc_w_cu];% go_plot_diag];
lg=legend(leg_gos,sprintf('Posterior\npredictive median'),...
    sprintf('C.I. w/o code\nuncertainty'),...
    sprintf('C.I. with code\nuncertainty'),...
    'Location','northeast');
% lg.Position(1:2)=[.623 .725];
flushLegend(lg,'northeast');
lg.Box='off';
% lgpos = lg.Position;
% lg.Position = lgpos + [-.004 -.002 -.004 -.002];

    
%%% Save
set(h,'Color','white');
% export_fig 	 -png -m3 -painters
% saveas(h,'FIG_cost_grid_pareto_with_code_uncert.png');
figstr = 'FIG_cost_grid_pareto_bands';
set(h,'PaperPositionMode','auto')
% print(h,figstr,'-depsc','-r600')

%% WTA Pareto bands
clc ; clearvars -except dpath ; close all ;

%%% Load the results
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-28_CTO_costgrid_size30';
% locstr= ['C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\'...
%     'stored_data\2020-04-26_CTO_costgrid_size30'];
load(locstr);

   
mean_y = results{1}.settings.mean_y ; std_y = results{1}.settings.std_y ;

% Collect Cost_lambdas, and posterior mean and sds for costs, defl, rot, as
% well as upper and lower .05 quantiles
m=length(results)-1; % Store number of target cost_lambdas
cost_lambda = zeros(m,1); % This will store cost_lambdas
cred_level = 90; % Set desired level for credible bands (in %)
alpha = (100-cred_level)/100; % Convert cred_level to alpha level
pmo = zeros(m,2); % This will store posterior mean output of emulator
pdo = zeros(m,2); % ``'' median output
pso = zeros(m,2); % ``'' appropriate multiple of standard deviations
plo = zeros(m,2); % ``'' lower (alpha/2) quantile
puo = zeros(m,2); % ``'' upper (alpha/2) quantile
for ii = 1:m % This loop populates the above arrays
    output_means = results{ii}.model_output.means .* std_y + mean_y;
    output_sds = sqrt(results{ii}.model_output.vars) .* std_y;
    pmo(ii,:) = mean(output_means);
    pdo(ii,:) = quantile(output_means,0.5);
    pso(ii,:) = norminv(1-alpha/2) * ...
        mean(output_sds);
    plo(ii,:) = quantile(output_means,alpha/2);
    puo(ii,:) = quantile(output_means,1-alpha/2);
    cost(ii) = results{ii}.settings.obs_y(1,end) * std_y(end) + mean_y(end);
end
% Now we break the arrays up each into 2 vectors, one for each output
post_cost_mean = pmo(:,2);
post_defl_mean = pmo(:,1);
post_cost_median = pdo(:,2);
post_defl_median = pdo(:,1);
post_cost_sd = pso(:,2);
post_defl_sd = pso(:,1);
post_cost_lq = plo(:,2);
post_cost_uq = puo(:,2);
post_defl_lq = plo(:,1);
post_defl_uq = puo(:,1);
% Get quantiles plus code uncertainty
post_cost_uq_cu = post_cost_uq + post_cost_sd;
post_cost_lq_cu = post_cost_lq - post_cost_sd;
post_defl_uq_cu = post_defl_uq + post_defl_sd;
post_defl_lq_cu = post_defl_lq - post_defl_sd;
% Get lims for the two sets of plots
ylimrat=1.01;
ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];

%%% Begin figures
% Set alphas for two types of uncertainty
alpha_wcu = 0.5;  %with code uncertainty
alpha_wocu= 0.15; %without
h=figure('rend','painters','pos',[10 (10) 360 240]);
x = linspace(ylim_cost(1),ylim_cost(2)); % x fills the cost domain
[subplts , pos] = tight_subplot(1,2,0.175,[0.15 0.02],[0.11 0.01]);

% Now begin plot 1/2
axes(subplts(1));
% Get main curve
pdefl = pchip(post_cost_median,post_defl_median,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(post_cost_median,post_defl_uq,x);
pdefllq = pchip(post_cost_median,post_defl_lq,x);

%     pdefluq_code_uncert = pchip(post_cost_median,post_defl_uq_cu,x);
%     pdefllq_code_uncert = pchip(post_cost_median,post_defl_lq_cu,x);
%     f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
%     set(f,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
%     hold on;
%     f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
%         fliplr(pdefluq)],'k');
%     unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
%         fliplr(pdefllq)],'k');

% unc_wo_cu = fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
% set(unc_wo_cu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
hold on;

pdefluq_code_uncert = pchip(post_cost_median,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(post_cost_median,post_defl_lq_cu,x);
% fwocu=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k',...
%     'FaceAlpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
% set(fwocu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
hold on;
% fwcu=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
%     fliplr(pdefluq)],'k',...
%     'FaceAlpha',alpha_wcu,'EdgeAlpha',alpha_wcu);
% unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
%     fliplr(pdefllq)],'k',...
%     'FaceAlpha',alpha_wcu,'EdgeAlpha',alpha_wcu);

median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
% plot(x,pdefluq,':k',...
%      x,pdefllq,':k');
xl2=xlabel('Target cost');
ylabel('Deflection');
xlim(ylim_cost);
ylim(ylim_defl);
% Add points
%     plot(post_cost_median,post_defl_median,'o');
%     plot(post_cost_median,pmo(:,1),'.');


figpos = get(h,'pos');

% Add NSGA-II results
nsga2_res = ...
    plot(results{end}.final_obj(:,2),results{end}.final_obj(:,1),'b*');


% Now add a legend.
ylim(ylim+[0 0.03]);
leg_gos = [nsga2_res median_line];% go_plot_diag];
lg=legend(leg_gos,sprintf('NSGA-II results'),...
    sprintf('Posterior\npredictive median'),...
    'Location','northeast');
% lg.Position(1:2)=[.623 .725];
flushLegend(lg,'northeast');
lg.Box='off';
% lgpos = lg.Position;
% lg.Position = lgpos + [-.004 -.002 -.004 -.002];

% Now add second plot. Just copy code above and zoom in.
axes(subplts(2));
unc_wo_cu = fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(unc_wo_cu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
hold on;

pdefluq_code_uncert = pchip(post_cost_median,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(post_cost_median,post_defl_lq_cu,x);
fwocu=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k',...
    'FaceAlpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
% set(fwocu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
hold on;
fwcu=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
    fliplr(pdefluq)],'k',...
    'FaceAlpha',alpha_wcu,'EdgeAlpha',alpha_wcu);
unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
    fliplr(pdefllq)],'k',...
    'FaceAlpha',alpha_wcu,'EdgeAlpha',alpha_wcu);

median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
% plot(x,pdefluq,':k',...
%      x,pdefllq,':k');
xl2=xlabel('Target cost');
ylabel('Deflection');
xlim_zoom = [106.1 106.2];
ylim_zoom = [.80235 .80257];
xlim(xlim_zoom);
ylim(ylim_zoom);
% Add NSGA-II results
nsga2_res = ...
    plot(results{end}.final_obj(:,2),results{end}.final_obj(:,1),'b*');

% Now add a legend.
% ylim(ylim_zoom+[0 0.0003]);
leg_gos = [fwocu fwcu];% go_plot_diag];
lg=legend(leg_gos,...
    sprintf([int2str(cred_level) '%% C.I. without\ncode uncertainty']),...
    sprintf([int2str(cred_level) '%% C.I. with\ncode uncertainty']),...
    'Location','northeast');
% lg.Position(1:2)=[.623 .725];
flushLegend(lg,'northeast');
lg.Box='off';
% lgpos = lg.Position;
% lg.Position = lgpos + [-.004 -.002 -.004 -.002];


%%% Save
set(h,'Color','white');
% export_fig 	 -png -m3 -painters
% saveas(h,'FIG_cost_grid_pareto_with_code_uncert.png');
figstr = sprintf(['FIG_cost_grid_pareto_bands']);
set(h,'PaperPositionMode','auto')
print(h,figstr,'-depsc','-r600')


%% WTA Histogram2 of highest density regions of posterior distribution 
clc ; clearvars -except dpath ; close all ;

%%% Load the calibration results
% clearvars -except dpath res ; 
locstr = [dpath,'stored_data\'...
    '2020-04-26_CTO_size30'];
load(locstr);
samps = res.theta1(res.settings.burn_in:end,:) ;

fighist = figure('pos',[10 10  360.0000  240]);
binwidth = 1.5e-3*[0.002 1];
h2=histogram2(samps(:,1),samps(:,2),...
    'Normalization','pdf','FaceAlpha',.75);
set(gca,'View',[225 30]);
set(gca,'ZTick',[]);
grid off;
ttlhist = title('Posterior distribution of \theta');
ttlhist.FontSize = 11;
xlbl = xlabel('Vol. fraction'); ylbl = ylabel('Thickness (mm)');
set(fighist,'Color','w');
% xlim([0.5 0.6]);
% ylim([10 10.6]);
% pause(1.5);

% Save it
figstr = 'FIG_post_dist_hist2';
set(fighist,'PaperPositionMode','auto')
ylbl.Units = 'pixels'; xlbl.Units='pixels';
ylbl.Position = ylbl.Position + [-8.0 2 0];
xlbl.Position = xlbl.Position + [8.0 2 0];
% print(fighist,figstr,'-depsc','-r600')


%% WTA prior predictive distribution vs posterior predictive distribution
clc ; clearvars -except dpath ; close all ;

%%% Load calib results
% clearvars -except dpath res ; close all;
% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO'];
% locstr = [dpath,'stored_data\'...
%     '2020-04-26_CTO_size30'];
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-30_CTO_size30';
load(locstr);
mean_y = res.settings.mean_y ; std_y = res.settings.std_y;
posamps = res.model_output.means .* std_y + mean_y;
des_obs = res.settings.obs_y(1,:).* std_y + mean_y;

%%% Load prior predictive results
locstr2 = [dpath,'stored_data\'...
    '2019-11-06_prior_predictive_distributions'];
load(locstr2);
prsamps = prior_model_output.means.* std_y + mean_y;


%%% Make figure using histograms
f=figure('pos',[10 10  360.0000  200]);
[subplts,pos] = tight_subplot(1,3,0.02,[ 0.08 0.01],0.03);
% Deflection
axes(subplts(1));
[p,x,bw]=ksdensity(posamps(:,1));
max_lim = max(p);
plot(x,p,'LineWidth',2);
set(gca,'YTick',[]);
%histogram(posamps(:,1),'Normalization','pdf','Edgecolor','none'); 
hold on;
[p,x]=ksdensity(prsamps(:,1));
max_lim = max([p(:);max_lim]);
plot(x,p,'--','LineWidth',2);
%histogram(prsamps(:,1),'Normalization','pdf','Edgecolor','none');
text(0.005,.9,...
    'Deflection','VerticalAlignment','bottom','Units','normalized');
text(1.715,102,'Rotation','VerticalAlignment','bottom');
% xlim([0.6 0.85]);
% ylim([0 110]);
line([des_obs(1) des_obs(1)],ylim,'Color','black','Linestyle',':',...
    'linewidth',2);
xlim([min(x)/1.075,max(x)]);

% Rotation
axes(subplts(2));
[p,x,bw]=ksdensity(posamps(:,2));
max_lim = max(p);
plot(x,p,'LineWidth',2);
set(gca,'YTick',[]);
%histogram(posamps(:,2),'Normalization','pdf','Edgecolor','none'); 
hold on;
[p,x]=ksdensity(prsamps(:,2));
max_lim = max([p(:);max_lim]);
plot(x,p,'--','LineWidth',2);
%histogram(prsamps(:,2),'Normalization','pdf','Edgecolor','none');
text(0.005,.9,...
    'Rotation','VerticalAlignment','bottom','Units','normalized');
% xlim([0.075,0.105])
% ylim([0 700]);
ylim = get(gca,'ylim');
line([des_obs(2) des_obs(2)],ylim,'Color','black','Linestyle',':',...
    'linewidth',2);
set(gca,'ylim',ylim);
xlim([min(x)/1.02,max(x)]);


% Cost
axes(subplts(3));
[p,x,bw]=ksdensity(posamps(:,3));
max_lim = max(p);
plot(x,p,'LineWidth',2);
set(gca,'YTick',[]);
%histogram(posamps(:,3),'Normalization','pdf','Edgecolor','none'); 
hold on;
[p,x]=ksdensity(prsamps(:,3));
max_lim = max([p(:);max_lim]);
plot(x,p,'--','LineWidth',2);
%histogram(prsamps(:,3),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Cost','VerticalAlignment','bottom','Units','normalized');
% ylim([0 .0700]);
xlim([60 400]);
ylim = [0 max_lim*1.15];%get(gca,'ylim');ylim
line([des_obs(3) des_obs(3)],ylim,'Color','black','Linestyle',':',...
    'linewidth',2);
set(gca,'ylim',ylim);


% Add suptitle
% st=suptitle('Prior and posterior predictive distributions');
% st.Position=[0.5 -.1 0];
[lg,icons,~,~]=legend('Posterior','Prior','Target','Location','northeast');
% flushLegend(lg,'northeast');
resizeLegend();
pos=lg.Position; 
lg.Position = pos + [.0915 0.055 0 0];
% 
% %%% Save
set(f, 'Color','white');
% % export_fig FIG_prior_vs_posterior_dist -eps -m3 -painters
figstr = 'FIG_prior_vs_posterior_dist';
set(f,'PaperPositionMode','auto')
print(f,figstr,'-depsc','-r600')


%% WTA_alt Pareto bands
clc ; clearvars -except dpath ; close all ;

% %%% Load the results
% locstr = [dpath,'stored_data\'...
%     '2020-04-23_CTO_costgrid_new_cases_size100'];
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-29_3_CTO_grids_size20';
load(locstr);

for hh=1:length(all_results)
    
    results=all_results{hh};
    
    mean_y = results{1}.settings.mean_y ; std_y = results{1}.settings.std_y ;

    % Collect Cost_lambdas, and posterior mean and sds for costs, defl, rot, as
    % well as upper and lower .05 quantiles
    m=length(results)-1; % Store number of target cost_lambdas
    cost_lambda = zeros(m,1); % This will store cost_lambdas
    cred_level = 90; % Set desired level for credible bands (in %)
    alpha = (100-cred_level)/100; % Convert cred_level to alpha level
    pmo = zeros(m,2); % This will store posterior mean output of emulator
    pdo = zeros(m,2); % ``'' median output
    pso = zeros(m,2); % ``'' appropriate multiple of standard deviations
    plo = zeros(m,2); % ``'' lower (alpha/2) quantile
    puo = zeros(m,2); % ``'' upper (alpha/2) quantile
    for ii = 1:m % This loop populates the above arrays
        burn_in = results{ii}.settings.burn_in +1;
        output_means = results{ii}.model_output.means(burn_in:end,:) .* std_y + mean_y;
        output_sds = sqrt(results{ii}.model_output.vars(burn_in:end,:)) .* std_y;
        pmo(ii,:) = mean(output_means);
        pdo(ii,:) = quantile(output_means,0.5);
        pso(ii,:) = norminv(1-alpha/2) * ...
            mean(output_sds);
        plo(ii,:) = quantile(output_means,alpha/2);
        puo(ii,:) = quantile(output_means,1-alpha/2);
        cost(ii) = results{ii}.settings.obs_y(1,end) * std_y(end) + mean_y(end);
    end
    % Now we break the arrays up each into 2 vectors, one for each output
    post_cost_mean = pmo(:,2);
    post_defl_mean = pmo(:,1);
    post_cost_median = pdo(:,2);
    post_defl_median = pdo(:,1);
    post_cost_sd = pso(:,2);
    post_defl_sd = pso(:,1);
    post_cost_lq = plo(:,2);
    post_cost_uq = puo(:,2);
    post_defl_lq = plo(:,1);
    post_defl_uq = puo(:,1);
    % Get quantiles plus code uncertainty
    post_cost_uq_cu = post_cost_uq + post_cost_sd;
    post_cost_lq_cu = post_cost_lq - post_cost_sd;
    post_defl_uq_cu = post_defl_uq + post_defl_sd;
    post_defl_lq_cu = post_defl_lq - post_defl_sd;
    % Get lims for the two sets of plots
    ylimrat=1.01;
    ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
    ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];

    %%% Begin figures
    % Set alphas for two types of uncertainty
    alpha_wcu = 0.5;  %with code uncertainty
    alpha_wocu= 0.15; %without
    h=figure('rend','painters','pos',[10 (10+(hh-1)*180) 360 110]);
    x = linspace(ylim_cost(1),ylim_cost(2)); % x fills the cost domain
%     [subplts , pos] = tight_subplot(1,2,0.175,[0.15 0.02],[0.11 0.01]);

    % Now begin plot 1/2
%     axes(subplts(1));
    % Get main curve
    pdefl = pchip(post_cost_median,post_defl_median,x);
    % Get upper and lower 0.05 quantiles curves
    pdefluq = pchip(post_cost_median,post_defl_uq,x);
    pdefllq = pchip(post_cost_median,post_defl_lq,x);
    
    %     pdefluq_code_uncert = pchip(post_cost_median,post_defl_uq_cu,x);
%     pdefllq_code_uncert = pchip(post_cost_median,post_defl_lq_cu,x);
%     f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
%     set(f,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
%     hold on;
%     f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
%         fliplr(pdefluq)],'k');
%     unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
%         fliplr(pdefllq)],'k');
    
    unc_wo_cu = fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
    set(unc_wo_cu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
    hold on;
    
    pdefluq_code_uncert = pchip(post_cost_median,post_defl_uq_cu,x);
    pdefllq_code_uncert = pchip(post_cost_median,post_defl_lq_cu,x);
    fwocu=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
    set(fwocu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
    hold on;
    fwcu=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
        fliplr(pdefluq)],'k');
    unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
        fliplr(pdefllq)],'k');
    
    median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
    % plot(x,pdefluq,':k',...
    %      x,pdefllq,':k');
    xl2=xlabel('Target cost');
    ylabel('Deflection');
    xlim(ylim_cost);
    ylim(ylim_defl);
    % Add points
%     plot(post_cost_median,post_defl_median,'o');
%     plot(post_cost_median,pmo(:,1),'.');
    

    figpos = get(h,'pos');

    % Add NSGA-II results
    nsga2_res = ...
        plot(results{end}.final_obj(:,2),results{end}.final_obj(:,1),'b*');


    % Now add a legend.
    ylim(ylim+[0 0.03]);
    leg_gos = [nsga2_res median_line];% go_plot_diag];
    lg=legend(leg_gos,sprintf('NSGA-II results'),...
        sprintf('Posterior\npredictive median'),...
        'Location','northeast');
    % lg.Position(1:2)=[.623 .725];
    flushLegend(lg,'northeast');
    lg.Box='off';
    % lgpos = lg.Position;
    % lg.Position = lgpos + [-.004 -.002 -.004 -.002];


    % saveas(h,'FIG_cost_grid_pareto.png');

    % Now add in code uncertainty. That is, the above assumes that the GP
    % emulator nails the FE code precisely. But of course the GP emulator has
    % nonnegligible variance. That's the code uncertainty. So our confidence
    % bands should reflect it. So we add it in here, by dropping the
    % appropriate multiple of the sd from each lower quantile and adding it to
    % each upper quantile.
    % First, open the figure prior to calling suptitle.
    % h=openfig('tempfig');
%     axes(subplts(2));

%     hold on;
%     median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
%     set(f,'facealpha',alpha_wcu,'EdgeAlpha',alpha_wcu);
%     set(unc_w_cu,'facealpha',alpha_wcu,'EdgeAlpha',alpha_wcu);
%     ylabel('Deflection');
%     xl2=xlabel('Target cost');
% %     xlim([183.7013 186.1534]);
% %     ylim([0.7255 0.7273]);
% 
%     % Add NSGA-II results
%     plot(results{end}.final_obj(:,2),results{end}.final_obj(:,1),'*');
% 
%     % Now add a main title and fix any infelicities
%     % suptitle(['Posterior estimate vs. target cost,',...
%     %     ' with ',num2str(cred_level),'% credible interval ']); 
%     set(h,'pos',figpos); % Just so we can reuse the positioning code from above
%     % p = get(xl1,'position');
%     % set(xl1,'position',p + [0 2.75 0]);
%     % p = get(xl2,'position');
%     % set(xl2,'position',p + [0 0.00125 0])
%     % p = get(xl3,'position');
%     % set(xl3,'position',p + [0 0.0002 0])
% 
%     % Now add a legend.
%     ylim(ylim+[0 0.0007]);
%     leg_gos = [median_line unc_wo_cu unc_w_cu];% go_plot_diag];
%     lg=legend(leg_gos,sprintf('Posterior\npredictive median'),...
%         sprintf('C.I. w/o code\nuncertainty'),...
%         sprintf('C.I. with code\nuncertainty'),...
%         'Location','northeast');
%     % lg.Position(1:2)=[.623 .725];
%     flushLegend(lg,'northeast');
%     lg.Box='off';
%     % lgpos = lg.Position;
%     % lg.Position = lgpos + [-.004 -.002 -.004 -.002];


    %%% Save
    set(h,'Color','white');
    % export_fig 	 -png -m3 -painters
    % saveas(h,'FIG_cost_grid_pareto_with_code_uncert.png');
    figstr = sprintf(['FIG_cost_grid_pareto_bands',int2str(hh)]);
    set(h,'PaperPositionMode','auto')
    print(h,figstr,'-depsc','-r600')

end

%% WTA_alt Pareto bands - double grid
clc ; clearvars -except dpath ; close all ;

% %%% Load the results
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-26_CTO_double_grid_size20';
load(locstr);

for hh=1:length(all_results)
    
    results=all_results{hh};
    
    mean_y = results{1}.settings.mean_y ; std_y = results{1}.settings.std_y ;

    % Collect Cost_lambdas, and posterior mean and sds for costs, defl, rot, as
    % well as upper and lower .05 quantiles
    m=length(results)-1; % Store number of target costs
    cred_level = 90; % Set desired level for credible bands (in %)
    alpha = (100-cred_level)/100; % Convert cred_level to alpha level
    pmo = zeros(m,2); % This will store posterior mean output of emulator
    pdo = zeros(m,2); % ``'' median output
    pso = zeros(m,2); % ``'' appropriate multiple of standard deviations
    plo = zeros(m,2); % ``'' lower (alpha/2) quantile
    puo = zeros(m,2); % ``'' upper (alpha/2) quantile
    for ii = 1:(m/2) % This loop populates the above arrays
        output_means = results{ii}.model_output.means .* std_y + mean_y;
        output_sds = sqrt(results{ii}.model_output.vars) .* std_y;
        pmo(ii,:) = mean(output_means);
        pdo(ii,:) = quantile(output_means,0.5);
        pso(ii,:) = norminv(1-alpha/2) * ...
            mean(output_sds);
        plo(ii,:) = quantile(output_means,alpha/2);
        puo(ii,:) = quantile(output_means,1-alpha/2);
        cost(ii) = results{ii}.settings.obs_y(1,end) * std_y(end) + mean_y(end);
    end
    for ii=(m/2+1):m
        output_means = results{ii}.model_output.means .* std_y + mean_y;
        output_sds = sqrt(results{ii}.model_output.vars) .* std_y;
        pmo(ii,:) = mean(output_means);
        pdo(ii,:) = quantile(output_means,0.5);
        pso(ii,:) = norminv(1-alpha/2) * ...
            mean(output_sds);
        plo(ii,:) = quantile(output_means,alpha/2);
        puo(ii,:) = quantile(output_means,1-alpha/2);
        defl(ii) = results{ii}.settings.obs_y(1,end) * std_y(end) + mean_y(end);
    end
    % Now we break the arrays up each into 2 vectors, one for each output
    post_cost_mean = pmo(:,2);
    post_defl_mean = pmo(:,1);
    post_cost_median = pdo(:,2);
    post_defl_median = pdo(:,1);
    post_cost_sd = pso(:,2);
    post_defl_sd = pso(:,1);
    post_cost_lq = plo(:,2);
    post_cost_uq = puo(:,2);
    post_defl_lq = plo(:,1);
    post_defl_uq = puo(:,1);
    % Get quantiles plus code uncertainty
    post_cost_uq_cu = post_cost_uq + post_cost_sd;
    post_cost_lq_cu = post_cost_lq - post_cost_sd;
    post_defl_uq_cu = post_defl_uq + post_defl_sd;
    post_defl_lq_cu = post_defl_lq - post_defl_sd;
    % Get lims for the two sets of plots
    ylimrat=1.01;
    ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
    ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];

    %%% Begin figures
    % Set alphas for two types of uncertainty
    alpha_wcu = 0.5;  %with code uncertainty
    alpha_wocu= 0.15; %without
    h=figure('rend','painters','pos',[10 (10+(hh-1)*180) 360 110]);
    x = linspace(ylim_cost(1),ylim_cost(2)); % x fills the cost domain
%     [subplts , pos] = tight_subplot(1,2,0.175,[0.15 0.02],[0.11 0.01]);

    % Now begin plot 1/2
%     axes(subplts(1));
    % Get main curve
    pdefl = pchip(post_cost_median,post_defl_median,x);
    % Get upper and lower 0.05 quantiles curves
    pdefluq = pchip(post_cost_median,post_defl_uq,x);
    pdefllq = pchip(post_cost_median,post_defl_lq,x);
        
    unc_wo_cu = fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
    set(unc_wo_cu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
    hold on;
    
    pdefluq_code_uncert = pchip(post_cost_median,post_defl_uq_cu,x);
    pdefllq_code_uncert = pchip(post_cost_median,post_defl_lq_cu,x);
    fwocu=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
    set(fwocu,'facealpha',alpha_wocu,'EdgeAlpha',alpha_wocu);
    hold on;
    fwcu=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
        fliplr(pdefluq)],'k');
    unc_w_cu=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
        fliplr(pdefllq)],'k');
    
    median_line = plot(x,pdefl,'-r','LineWidth',1.5); % Mean
    % plot(x,pdefluq,':k',...
    %      x,pdefllq,':k');
    xl2=xlabel('Target cost');
    ylabel('Deflection');
    xlim(ylim_cost);
    ylim(ylim_defl);
    % Add points
    plot(post_cost_median,post_defl_median,'o');
    plot(post_cost_median,pmo(:,1),'.');
    

    figpos = get(h,'pos');

    % Add NSGA-II results
    nsga2_res = ...
        plot(results{end}.final_obj(:,2),results{end}.final_obj(:,1),'b*');


    % Now add a legend.
    ylim(ylim+[0 0.03]);
    leg_gos = [nsga2_res median_line];% go_plot_diag];
    lg=legend(leg_gos,sprintf('NSGA-II results'),...
        sprintf('Posterior\npredictive median'),...
        'Location','northeast');
    % lg.Position(1:2)=[.623 .725];
    flushLegend(lg,'northeast');
    lg.Box='off';
    % lgpos = lg.Position;
    % lg.Position = lgpos + [-.004 -.002 -.004 -.002];


    %%% Save
    set(h,'Color','white');
    % export_fig 	 -png -m3 -painters
    % saveas(h,'FIG_cost_grid_pareto_with_code_uncert.png');
    figstr = sprintf(['FIG_cost_grid_pareto_bands',int2str(hh)]);
    set(h,'PaperPositionMode','auto')
    print(h,figstr,'-depsc','-r600')

end


%% ZDT Figure describing the problem
clc ; clearvars -except dpath ; close all; 
% Take a look at surfaces of example function output
theta1=linspace(0,1);
theta2=linspace(0,1);
% We need to convert these to a two-col matrix of all possible combinations
[ T1 , T2 ] = meshgrid(theta1,theta2); 
ths = cat(2,T1',T2');
Theta = reshape(ths,[],2);

% Set other parameters
thetas = zeros(size(Theta,1),3);

% Now get output values
opc = TP_ZDT1_objfun([Theta thetas]);

% Normalize the outputs
oscl_n = (opc(:,1)-min(opc(:,1))) / range(opc(:,1)) ;
perf_n = (opc(:,2)-min(opc(:,2))) / range(opc(:,2)) ;

% Now take a look at the surfaces
oscls = reshape(oscl_n,[],length(theta1));
perfs = reshape(perf_n,[],length(theta1));

e1 = [  0   0   0 ] ;  % output 1 color
e2 = [.4 .4 .4 ] ;  % output 2 color

g1 = [ 0 0 .6 ] ; 
g2 = [ 1 .3 0 ] ;
g3 = [ .4 1 .8 ] ; 
f1 = 'r'; f2='g'; f3='b'; % for color version
ea = 1      ;  % edge alpha
fa = 0      ;  % face alpha


f=figure('pos',[10 10 325 200]);
sp1 = subplot(1,2,1);
sp2 = subplot(1,2,2);
set(sp1,'position',[.1 .225 .65 .825]);
set(sp2,'position',[.4 .225 .6 .825]);
set(f,'currentaxes',sp1);
surf(theta2,theta1,oscls,'FaceColor',g1,'EdgeColor',g1,...
    'EdgeAlpha',ea,'FaceAlpha',fa);
% axis vis3d;

hold on;
surf(theta2,theta1,perfs,'FaceColor',g2,'EdgeColor',g2,...
    'EdgeAlpha',ea,'FaceAlpha',fa); 
% axis vis3d;

% axis vis3d;
xlabel('\theta_2'); ylabel('\theta_1'); zlabel('Outcomes');

h=gca;
h.View = [97 15] ; % Sets the perspective
h.View = [108 4];
set(h,'ztick',[]);

% title('Model outcomes on normalized scale');

% hh = legend('y_1','y_2','y_3','Orientation','horizontal',...
%     'Location','south');
% hh.Position = hh.Position + [ 0 -.115 0 0 ];
% set(gcf,'unit','inches');
% fig_size = get(gcf,'position');
set(f,'currentaxes',sp2);
s1=patch(nan*[0 1 1 0],[0 0 1 1],g1);
hold on;
s2=patch(nan*[1 2 2 1],[1 1 2 2],g2);
s3=patch(nan*[2 3 3 2],[2 2 3 3],g3);
set(sp2,'Visible','off');
% set(s1,'Visible','off');set(s2,'Visible','off');set(s3,'Visible','off');

hh = legend(sp2,'y_1','y_2','Orientation','vertical',...
    'Location','east');
% set(hh,'unit','inches');
% leg_size = get(hh,'position');
% fig_size(3) = fig_size(3) + leg_size(3);
% set(gcf,'position',fig_size)

% %%% Code to turn it into the version used on SCSC poster
% hh.Orientation = 'vertical';
% hh.Location = 'east';
% f.Position = [360.3333  197.6667  452.0000  314.6667];
% set(f,'Color',[251/255 244/255 245/255]);
% set(h,'Color',[251/255 244/255 245/255]);
% set(hh,'Color',[251/255 244/255 245/255]);

%saveas(f,'FIG_toy_sim_model_outputs.png');
figstr = 'FIG_ZDT_model_outputs';
set(f,'Color','w');
% export_fig FIG_toy_sim_model_outputs -eps -m3 -painters
print(f,figstr,'-depsc','-r600');

%% ZDT Results (output on contour plot) with arbitrary target

clc ; clearvars -except dpath ; close all ;

% Initialize figure
figposition = [10 10 390 260];
fig1 = figure('Position',figposition);
ttlfontsize = 11;

% Load results
% locstr = [dpath,'Example\Ex_results\'...
%     '2020_04_19_CTO_noemulator_TP_ZDT1_5p222SD'];
% locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-27_TP_ZDT1_5p222SD';
% locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-28_TP_ZDT1_5p222SD';
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-29_TP_ZDT1_5p222SD';
load(locstr);

% Get arbitrary desired observation
mean_y = res.settings.mean_y;
std_y   = res.settings.std_y;
des_obs = res.settings.obs_y; 

% Set true optimum: this is for the arbitrary desired observation specified
% above. Found through fmincon (elsewhere).
optim = [ 0.4958 0 0 0 0] ;

burn_in = res.settings.burn_in;
samps = [res.theta1(burn_in:end,1) res.theta1(burn_in:end,2)];
[theta1,theta2] = meshgrid(linspace(0,1,1000),linspace(0,1,1000));
% Get true samples from a dense mesh grid over the supports of t1 and t2 at
% optimal values of t3-5;
M=1000;
theta1=linspace(0,1,M);
theta2=linspace(0,1,M);
[t1, t2] = meshgrid(theta1,theta2);
t2_5val = 0;
all_ctheta = [ t1(:) t2(:) t2_5val*ones(M^2,1) t2_5val*ones(M^2,1) t2_5val*ones(M^2,1)] ;
output_over_mesh = TP_ZDT1_objfun(all_ctheta);

%%% Put the outputs and the desired observation on the standardized scale
outputs_std = (output_over_mesh - mean_y)./std_y;

% Now get Euclidean norms of each standardized output
dists = sqrt ( sum ( (outputs_std-des_obs).^2 , 2 ) ) ;
redists = reshape(dists,1000,1000);

% Get distance of optim output from des_obs
output_at_optim_std = (TP_ZDT1_objfun(optim)-mean_y)./std_y;
dist_optim = sqrt ( sum ( (output_at_optim_std-des_obs).^2 , 2 ) )

%%% Make scatterhist of posterior samples
colormap autumn
subplot(1,2,1);
sc1=scatterhist(samps(:,1),samps(:,2),'Marker','.','Color','b',...
    'Markersize',1,'Kernel','on'); 
% sc(2).Children.EdgeColor = 'none';
% sc(2).Children.FaceAlpha = 1;
% sc(3).Children.EdgeColor = 'none';
% sc(3).Children.FaceAlpha = 1;
hold on;
ttl1 = title({'Posterior \theta samples:' ...
    'target [0.25, -6]'});
ttl1.FontSize = ttlfontsize;

%%% Now add contour plot 
contour_arr = [5.25 5.5 5.75 6 6.25 6.5];
[C1,h1]= contour(theta1,theta2,redists,contour_arr,'LineWidth',3);
% [C,h]= contour(theta1,theta2,redists,[16 17 18 19 20 ],'LineWidth',3);
clabel(C1,h1,'fontsize',12);
xlabel('\theta_1'); ylabel('\theta_2');
xlim([0 1]) ; ylim([0 1]);

%%% Add true optimum
p=plot(optim(1),optim(2),'ok','MarkerSize',7,'MarkerFaceColor','m',...
    'LineWidth',2);




% Start second figure
fig2 = figure('Position',figposition);

% Load results
% locstr = [dpath,'Example\Ex_results\'...
%     '2020_04_19_CTO_noemulator_TP_ZDT1_2SD'];
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-29_TP_ZDT1_2SD';
load(locstr);

% Set desired observation as the utopia point of the system
des_obs_original = des_obs; % Store the original first
des_obs = res.settings.obs_y;

% Set true optimum: this is for the arbitrary desired observation specified
% above. Found through fmincon (elsewhere).
optim = [ 0.4958 0 0 0 0] ;

burn_in = res.settings.burn_in;
samps = [res.theta1(burn_in:end,1) res.theta1(burn_in:end,2)];
[theta1,theta2] = meshgrid(linspace(0,1,1000),linspace(0,1,1000));

% Now get Euclidean norms of each standardized output
dists = sqrt ( sum ( (outputs_std-des_obs_original).^2 , 2 ) ) ;
redists = reshape(dists,1000,1000);

% Get distance of optim output from des_obs
output_at_optim_std = (TP_ZDT1_objfun(optim)-mean_y)./std_y;
dist_optim = sqrt ( sum ( (output_at_optim_std-des_obs_original).^2 , 2 ) )

%%% Make scatterhist of posterior samples
colormap autumn
sc=scatterhist(samps(:,1),samps(:,2),'Marker','.','Color','b',...
    'Markersize',1,'Kernel','on'); 
% sc(2).Children.EdgeColor = 'none';
% sc(2).Children.FaceAlpha = 1;
% sc(3).Children.EdgeColor = 'none';
% sc(3).Children.FaceAlpha = 1;
hold on;
ttl2 = title({'Posterior \theta samples:', ...
    'target [-0.40, -2.15]'});
ttl2.FontSize = ttlfontsize;

%%% Now add contour plot 
contour_arr = [2 2.25 2.5 2.75 3 3.25 5.25 5.5 5.75 6 6.25 6.5];
[C,h]= contour(theta1,theta2,redists,contour_arr,'LineWidth',3);
% [C,h]= contour(theta1,theta2,redists,[16 17 18 19 20 ],'LineWidth',3);
clabel(C,h,'fontsize',12);
xlabel('\theta_1'); ylabel('\theta_2');
xlim([0 1]) ; ylim([0 1]);

%%% Add true optimum
p=plot(optim(1),optim(2),'ok','MarkerSize',7,'MarkerFaceColor','m',...
    'LineWidth',2);

% pause
pause(.5);

% Mess with color and position
fig1.Children(2).Children(1).Visible='off';
fig1.Children(3).Children(1).Visible='off';
fig2.Children(2).Children(1).Visible='off';
fig2.Children(3).Children(1).Visible='off';
set(fig1,'Color','white');
sppos1 = fig1.Children(4).Position;
fig1.Children(4).Position = sppos1 + [-0.025 -0.035 .025 0];
set(fig2,'Color','white');
sppos2 = fig2.Children(4).Position;
fig2.Children(4).Position = sppos2 + [-0.025 -0.035 .025 0];
pause(1); % Seems to be necessary
xpdfpos = fig1.Children(2).Position;
fig1.Children(2).Position = xpdfpos + [0 0 0 .015];
xpdfpos = fig2.Children(2).Position;
fig2.Children(2).Position = xpdfpos + [0 0 0 .015];


% Save figures as .eps
figstr1 = sprintf('FIG_ZDT_results_5p222SD');
figstr2 = sprintf('FIG_ZDT_results_2SD');
% export_fig(figstr1,'-eps','-q0','-painters',fig1);
% export_fig(figstr2,'-eps','-q0','-painters',fig2);
% print(fig1,figstr1,'-depsc','-r600');
% print(fig2,figstr2,'-depsc','-r600');


