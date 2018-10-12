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

%% Heatmaps for posteriors in 3-part iter. calib. on toy sim. example
clc ; clearvars -except dpath ; close all ; 

% Load results
load([dpath,'Example\Ex_results\'...
    '2018-09-12_iterative_calibration'],...
    'all_results');

des_obs = all_results{1}.settings.desired_obs;

h1 = calib_heatmap(des_obs, all_results{1}.samples_os);
pos = h1.Position; h1.Position = pos + [ 0 0 -1/4 * pos(3) 0 ]
pos = h1.Children(4).Position; dist = pos(4); sc=1/7;
set(h1.Children(4),'Position', pos + [-dist*sc -dist*sc dist*sc dist*sc]);
title('Calibration iteration 1/3, \lambda_\delta = 1/256');
xlabel('\theta_1');ylabel('\theta_2');

h2 = calib_heatmap(des_obs, all_results{2}.samples_os);
pos = h2.Position; h2.Position = pos + [ 0 0 -1/4 * pos(3) 0 ]
pos = h2.Children(4).Position; dist = pos(4); sc=1/7;
set(h2.Children(4),'Position', pos + [-dist*sc -dist*sc dist*sc dist*sc]);
title('Calibration iteration 2/3, \lambda_\delta = 16.7');
xlabel('\theta_1');ylabel('\theta_2');

h3 = calib_heatmap(des_obs, all_results{3}.samples_os);
pos = h3.Position; h3.Position = pos + [ 0 0 -1/4 * pos(3) 0 ]
pos = h3.Children(4).Position; dist = pos(4); sc=1/7;
set(h3.Children(4),'Position', pos + [-dist*sc -dist*sc dist*sc dist*sc]);
title('Calibration iteration 3/3, \lambda_\delta = 1784.5');
xlabel('\theta_1');ylabel('\theta_2');

% Save them
set(h1,'Color','white');
set(h2,'Color','white');
set(h3,'Color','white');
% export_fig('FIG_iter_calib_1o3','-png','-m3','-painters',h1);
% export_fig('FIG_iter_calib_2o3','-png','-m3','-painters',h2);
% export_fig('FIG_iter_calib_3o3','-png','-m3','-painters',h3);

%% Get marginal posteriors from toy sim. example at all three levels
clc ; clearvars -except dpath ; close all ;

% Load results
load([dpath,'Example\Ex_results\'...
    '2018-09-12_iterative_calibration'],...
    'all_results');
samps1 =all_results{1}.samples_os(all_results{1}.settings.burn_in+2:end,:);
samps2 =all_results{2}.samples_os(all_results{2}.settings.burn_in+2:end,:);
samps3 =all_results{3}.samples_os(all_results{3}.settings.burn_in+2:end,:);
samps=cell(3,1); samps{1}=samps1; samps{2}=samps2; samps{3}=samps3;

% Make figures
h1 = figure('rend','painters','pos',[10 10 610 190]) ;
h2 = figure('rend','painters','pos',[10 10 610 190]) ;
h3 = figure('rend','painters','pos',[10 10 610 190]) ;
h=cell(3,1); h{1}=h1; h{2}=h2; h{3}=h3;

% Load true optimum
load([dpath,'Example\Ex_results\'...
    '2018-09-11_true_optimum_wrt_0']);

for ii = 1 : 3

    %%% Get the marginal plots
    %h1 = figure('rend','painters','pos',[10 10 610 160]) ; 
    set (groot, 'currentfigure', h{ii});
    subplot(1,2,1);
    histogram(samps{ii}(:,1), 'Normalization','pdf') ;
    xlim([0 3]);
    unifval = 1/3;
    hold on;
    %xlabel('\theta_1');
    ylims=[0 3.3];%ylim;
    plot([0 3], [unifval unifval],'--r','LineWidth',1);
    plot([optim(1) optim(1)], [ylims],':g','LineWidth',1.5);
    ylim(ylims);
%     set(gca, 'XTick', sort([optim(1), get(gca, 'XTick')]));
    

    subplot(1,2,2);
    histogram(samps{ii}(:,2), 'Normalization','pdf') ;
    xlim([0 6]);
    unifval = 1/6;
    hold on;
    %xlabel('\theta_2');
    ylims=[0 12];%ylim;
    plot([0 6], [unifval unifval],'--r','LineWidth',1);
    plot([optim(2) optim(2)], [ylims],':g','LineWidth',1.5);
    ylim(ylims);
    

%     %%% Save
    set(h{ii},'Color','white');
    figstr = sprintf('FIG_iter_post_marginals_%go3',ii);
    export_fig(figstr,'-png','-m3','-painters',h{ii});
    
end

%% Marginal posteriors from wind turbine appl. at all three levels
clc ; clearvars -except dpath ; close all ;

load([dpath,'stored_data\'...
    '2018-09-12_iterative_calibration'],...
    'all_results');
samps1 =all_results{1}.samples_os(all_results{1}.settings.burn_in+2:end,:);
samps2 =all_results{2}.samples_os(all_results{2}.settings.burn_in+2:end,:);
samps3 =all_results{3}.samples_os(all_results{3}.settings.burn_in+2:end,:);
samps=cell(3,1); samps{1}=samps1; samps{2}=samps2; samps{3}=samps3;

% Make figures
h1 = figure('rend','painters','pos',[10 10 610 190]) ;
h2 = figure('rend','painters','pos',[10 10 610 190]) ;
h3 = figure('rend','painters','pos',[10 10 610 190]) ;
h=cell(3,1); h{1}=h1; h{2}=h2; h{3}=h3;

for ii = 1 : 3

    %%% Get the marginal plots
    %h1 = figure('rend','painters','pos',[10 10 610 160]) ; 
    set (groot, 'currentfigure', h{ii});
    subplot(1,2,1);
    histogram(samps{ii}(:,1), 'Normalization','pdf') ;
    xlim([0.2 0.6]);
    unifval = 1/.4;
    hold on;
    %xlabel('Volume fraction');
    ylims=[0 44];%ylim;
    plot([0.2 0.6], [unifval unifval],'--r','LineWidth',2);
    ylim(ylims);
%     set(gca, 'XTick', sort([optim(1), get(gca, 'XTick')]));
    

    subplot(1,2,2);
    histogram(samps{ii}(:,2), 'Normalization','pdf') ;
    xlim([10 25]);
    unifval = 1/15;
    hold on;
    %xlabel('Thickness (mm)');
    ylims=[0 1.16];%ylim;
    plot([10 25], [unifval unifval],'--r','LineWidth',2);
    ylim(ylims);
    

    %%% Save
    set(h{ii},'Color','white');
    figstr = sprintf('FIG_wta_iter_post_marginals_%go3',ii);
    export_fig(figstr,'-png','-m3','-painters',h{ii});
end

%% Marginal posterior predictive dists from all levels of WTA
clc ; clearvars -except dpath ; close all ;

clc ; clearvars -except dpath ; close all ;

%%% Load prior predictive results and all levels of calib
load([dpath,'stored_data\'...
    '2018-09-03_prior_pred_distrib'],...
    'prior_pred_dist');
prsamps = prior_pred_dist.prior_pred_pts;
clear prior_pred_dist;
load([dpath,'stored_data\'...
    '2018-09-12_iterative_calibration'],...
    'all_results');
posamps = cell(3,1);
posamps{1} = all_results{1}.model_output.by_sample_est(...
    all_results{1}.settings.burn_in+2:end,:) ;
posamps{2} = all_results{2}.model_output.by_sample_est(...
    all_results{2}.settings.burn_in+2:end,:) ;
posamps{3} = all_results{3}.model_output.by_sample_est(...
    all_results{3}.settings.burn_in+2:end,:) ;
des_obs = all_results{1}.settings.desired_obs;
des_obs_upd = all_results{3}.settings.desired_obs;
clear all_results;
f=cell(3,1);

%%% Loop through levels, creating a figure for each
for ii = 1:3
    
    f{ii}=figure('pos',[10 10  720.0000  180]);
    % Deflection
    subplot(1,3,1);
    histogram(posamps{ii}(:,1),'Normalization','pdf','Edgecolor','none'); 
    hold on;
    histogram(prsamps(:,1),'Normalization','pdf','Edgecolor','none');
    text(0.485,.875,...
        'Deflection','VerticalAlignment','bottom','Units','normalized');
%     text(1.715,102,'Rotation','VerticalAlignment','bottom');
%     title('Deflection');
    xlim([0.6 0.85]);
    ylim([0 225]);
    line([des_obs(1) des_obs(1)],ylim,'Color','black','Linestyle',':',...
        'linewidth',1);
    line([des_obs_upd(1) des_obs_upd(1)],ylim,...
        'Color','red','Linestyle','--','linewidth',1);
    % Rotation
    subplot(1,3,2);
    histogram(posamps{ii}(:,2),'Normalization','pdf','Edgecolor','none'); 
    hold on;
    histogram(prsamps(:,2),'Normalization','pdf','Edgecolor','none');
    text(0.55,.875,...
        'Rotation','VerticalAlignment','bottom','Units','normalized');
%     title('Rotation');
    xlim([0.075,0.105])
    ylim([0 1625]);
    line([des_obs(2) des_obs(2)],ylim,'Color','black','Linestyle',':',...
        'linewidth',1);
    line([des_obs_upd(2) des_obs_upd(2)],ylim,...
        'Color','red','Linestyle','--','linewidth',1);  
    % Cost
    subplot(1,3,3);
    histogram(posamps{ii}(:,3),'Normalization','pdf','Edgecolor','none'); 
    hold on;
    histogram(prsamps(:,3),'Normalization','pdf','Edgecolor','none');
    text(0.72,.875,...
        'Cost','VerticalAlignment','bottom','Units','normalized');
%     title('Cost');
    ylim([0 .1100]);
    line([des_obs(3) des_obs(3)],ylim,'Color','black','Linestyle',':',...
        'linewidth',1);
    line([des_obs_upd(3) des_obs_upd(3)],ylim,...
        'Color','red','Linestyle','--','linewidth',1);
    % Add suptitle
    %st=suptitle('Prior (red) and posterior (blue) predictive distributions');
    %st.Position=[0.5 -.1 0];

    %%% Save
    set(f{ii}, 'Color','white');
    figstr = sprintf('FIG_iter_prior_vs_posterior_dist_%go3',ii);
    export_fig(figstr,'-png', '-m3', '-painters');
    
end
