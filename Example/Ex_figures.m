% Simulation example figures
% Here is code for figures based on the toy simulation example for the
% calibration to desired observations


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

%% Get deflection, rotation "bands" as function of cost, w\ true trade-off
clc ; clearvars -except dpath ; close all;
% First using the one big calibration including cost
load([dpath,'Example\Ex_results\'...
'2018-05-17_d0_incl_min_cost'],...
'results'); % This is the results
load([dpath,'Example\Ex_results\'...
    '2018-05-18_full_calib_emulator_output_estimates'],...
    'emout'); % This is the GP estimate of output at each draw of the MCMC


% user defined values
alpha = 0.1;
cred_level = (1-alpha)*100;
k=1200;
x=linspace(15,30); % cost values at which to plot

% Get quantiles of MCMC samples at each whole dollar amount
costs = 15:30;
data_at_cost = cell(length(costs),1);
dir_data_at_cost = cell(length(costs),1);
medians = zeros(length(costs),2);
dir_medians = zeros(length(costs),2);
uqs = zeros(length(costs),2);
lqs = zeros(length(costs),2);
dir_uqs = zeros(length(costs),2);
dir_lqs = zeros(length(costs),2);
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');
data_true = ctheta_output_nondom(:,4:6);
data = emout.output_means;

for ii = 1:16
    
    data_at_cost{ii} = data(round(data(:,3))==costs(ii),1:2);
    medians(ii,:) = median(data_at_cost{ii});
    uqs(ii,:)     = quantile(data_at_cost{ii},1-alpha/2);
    lqs(ii,:)     = quantile(data_at_cost{ii},alpha/2);
    
    dir_data_at_cost{ii} = data_true(round(data_true(:,3))==costs(ii),1:2);
    dir_medians(ii,:) = median(dir_data_at_cost{ii});
    dir_uqs(ii,:)     = quantile(dir_data_at_cost{ii},1-alpha/2);
    dir_lqs(ii,:)     = quantile(dir_data_at_cost{ii},alpha/2);
    
end

% Plot "pareto bands" for defl, rot as fn of cost, each individually
% Get main curve
pdefl   = pchip(costs,medians(:,1),x);
pdefluq = pchip(costs,uqs(:,1),x);
pdefllq = pchip(costs,lqs(:,1),x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
hold on;
set(f,'facealpha',.25);
plot(x,pdefl,'r'); 
plot(x,pdefluq,'k'); plot(x,pdefllq,'k');
xlabel('Cost'); ylabel('Deflection');
title(['Deflection vs. target cost,',...
    ' with ',num2str(cred_level),'% credible interval']); 
% Add true vals
% Modify stuff since nondoms stop at $22
costs = 15:21;
x = linspace(15,21);
dir_pdefl   = pchip(costs,dir_medians(1:7,1),x);
dir_pdefluq = pchip(costs,dir_uqs(1:7,1),x);
dir_pdefllq = pchip(costs,dir_lqs(1:7,1),x);
f=fill([ x , fliplr(x) ], [dir_pdefluq, fliplr(dir_pdefllq)],'g');
hold on;
set(f,'facealpha',.25);
plot(x,dir_pdefl,'b'); 
plot(x,dir_pdefluq,'g'); plot(x,dir_pdefllq,'g');
xlabel('Cost'); ylabel('Deflection');
title(['Deflection vs. target cost,',...
    ' with ',num2str(cred_level),'% credible interval']); 


% Now rotation
h=figure('rend','painters','pos',[10 10 500 320]);
protn   = pchip(costs,medians(:,2),x);
protnuq = pchip(costs,uqs(:,2),x);
protnlq = pchip(costs,lqs(:,2),x);
f=fill([ x , fliplr(x) ], [protnuq, fliplr(protnlq)],'k');
hold on;
set(f,'facealpha',.25);
plot(x,protn,'r'); 
plot(x,protnuq,'k'); plot(x,protnlq,'k');
xlabel('Cost'); ylabel('Rotation');
title(['Rotation vs. target cost,',...
    ' with ',num2str(cred_level),'% credible interval']); 



% Get median of MCMC samples and plot
median_vals = nearquant(.5,x,emout.output_means(:,1),...
    emout.output_means(:,3),k);
plot(x,median_vals,'-o');

% Get median of direct data and plot
% ( Get direct data from Ex_optimization_workspace, allperfs )
load([dpath,'Example\Ex_results\'...
    '2018-05-18_direct_data_nondom'],...
    'nondom_ap');

%% Get figure at particular cost describing rotation/deflection tradeoff
clc; clearvars -except dpath ; close all;
% load([dpath,'Example\Ex_results\'...
% '2018-05-17_d0_incl_min_cost'],...
% 'results'); % This is the results
load([dpath,'Example\Ex_results\'...
'2018-05-29_set_obs_var_d0'],...
'results'); % This is the results

% load([dpath,'Example\Ex_results\'...
%     '2018-05-18_full_calib_emulator_output_estimates'],...
%     'emout'); % This is the GP est of output at each draw of the MCMC

cost_target = 20;

% Get outputs at that cost
% pbi_samps = emout.output_means(results.settings.burn_in:end,:);
pbi_samps = ...
    results.model_output.by_sample_true(results.settings.burn_in:end,:);
samps_at_cost_idx = ...
    (round(pbi_samps(:,3)) ...
    == cost_target);
samps_at_cost = ...
    pbi_samps(samps_at_cost_idx,:);

% Take a look
plot(samps_at_cost(:,1),samps_at_cost(:,2),'ob');

% Now find the same using direct data
% load([dpath,'Example\Ex_results\'...
%     '2018-05-25_true_ctheta-output'],...
%     'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output'],...
    'ctheta_output');
dd_at_cost_idx = (round(ctheta_output(:,6)) == cost_target);
dd_at_cost = ctheta_output(dd_at_cost_idx,4:5);

% Take a look at direct data tradeoff
hold on;
plot(dd_at_cost(:,1),dd_at_cost(:,2),'og');

% Look specifically at nondominated ones
[nd_dd_at_cost, nd_dd_at_cost_idx] = nondominated(dd_at_cost);
plot(nd_dd_at_cost(:,1),nd_dd_at_cost(:,2),'or');


%% Get figure describing the problem
% Take a look at surfaces of example function output
theta1=linspace(0,3);
theta2=linspace(0,6);
% We need to convert these to a two-col matrix of all possible combinations
[ T1 , T2 ] = meshgrid(theta1,theta2); 
ths = cat(2,T1',T2');
Theta = reshape(ths,[],2);

% Set other parameters
c = repmat(2.5,size(Theta,1),1);

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

%% Get 3d figure of nondominated direct data plus nondom full calib MCMC
% Both outputs and inputs. In that order.
clc; clearvars -except dpath ; close all;
% load([dpath,'Example\Ex_results\'...
% '2018-05-17_d0_incl_min_cost'],...
% 'results');
% Use new results from set total observation variance method:
load([dpath,'Example\Ex_results\'...
'2018-05-28_d0_incl_min_cost'],...
'results');
samps_os = results.samples_os(results.settings.burn_in+2:end,:);

% load([dpath,'Example\Ex_results\'...
%     '2018-05-18_full_calib_emulator_output_estimates'],...
%     'emout');
emout.output_means = results.model_output.by_sample_est(...
    results.settings.burn_in+2:end,:);

% Use GP estimates rather than true function:
[nondoms,ndidx] = nondominated(emout.output_means);
nondom_inputs = samps_os(ndidx,:);
% ( Get direct data from Ex_optimization_workspace, allperfs )
%nondom_ap = nondominated(allperfs);
% Old version of Ex_sim direct data:
% load([dpath,'Example\Ex_results\'...
%     '2018-05-25_true_ctheta-output_nondom'],...
%     'ctheta_output_nondom');
% Current version:
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');
nondom_ap = ctheta_output_nondom(:,4:6);
nondom_ddin = ctheta_output_nondom(:,2:3);

% Take a look
Circlesize=20;
figure; h1 = scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    1,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',0.3,...
    'MarkerEdgeAlpha',0.3);
axis vis3d;
hold on;
scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0.4);
xlabel('Cost'); ylabel('Deflection'); zlabel('Rotation');
title(['Nondominated direct data (red) and MCMC sample est. output (blue)'...
    '']);
viewpt = [45 20];
view(viewpt);
gif('FIG_nd.gif','frame',gcf);
nfms = 120;
for ii = 1:nfms
    viewpt = viewpt + [ 360/nfms 0 ];
    view(viewpt);
    gif
end

% Now make the same figure, but look at the true output at mcmc sample pts
% Use the true function to find the output at these sample points
y_samps_true = Ex_sim( [2*ones(size(samps_os,1),1) samps_os]);
% Get nondominated outcomes
[nondoms,ndidx] = nondominated(y_samps_true);
figure; h1 = scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    1,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',0.3,...
    'MarkerEdgeAlpha',0.3);
axis vis3d;
hold on;
scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0.4);
xlabel('Cost'); ylabel('Deflection'); zlabel('Rotation');
title(['Nondominated direct data (red) and MCMC sample true output (blue)'...
    '']);
viewpt = [45 20];
view(viewpt);
gif('FIG_nd_true.gif','frame',gcf);
nfms = 120;
for ii = 1:nfms
    viewpt = viewpt + [ 360/nfms 0 ];
    view(viewpt);
    gif
end

% Now get a comparison of inputs. First, put up the nondominated inputs:
figure; h2 = scatter(nondom_ddin(:,1),nondom_ddin(:,2),Circlesize,'r',...
    'filled');
hold on;
scatter(nondom_inputs(:,1),nondom_inputs(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',0.25);
xlim([0 3]); ylim([0 6]);
title('MCMC samples filtered by estimated nondominance');
saveas(h2,'FIG_nd_in');

%% Get 3d figure of nondominated direct data plus nondom cost grid MCMC
clc ; clearvars -except dpath ; close all;
load([dpath,'Example\Ex_results\'...
'2018-05-17_cost_grid_results'],...
'results');

emout = []; % This will hold all outputs from all points in cost grid
for ii = 1:size(results,1)
    emout = [ emout ; results{ii}.model_output.by_sample ];
end

% Use GP estimates rather than true function to find nondominated samps:
nondoms = nondominated(emout);
% ( Get direct data from Ex_optimization_workspace, allperfs )
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');
nondom_ap = ctheta_output_nondom(:,4:6);

% Take a look
Circlesize=50;
%figure; h1 = scatter3(y_samps(:,3),y_samps(:,1),y_samps(:,2),...
%    Circlesize,'b','filled','MarkerFaceAlpha',.4);
figure; h1 = scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    Circlesize,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',.2);
axis vis3d;
hold on;
scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',.4);
xlabel('Cost'); ylabel('Deflection'); zlabel('Rotation');
title(['Nondominated direct data (red) and MCMC samples'...
    '(blue) from cost grid']);

%% 3d figure of nondom'd direct data w/ all samps output from full calib
clc; clearvars -except dpath ; close all;
load([dpath,'Example\Ex_results\'...
'2018-05-17_d0_incl_min_cost'],...
'results');
samps_os = results.samples_os(results.settings.burn_in:end,:);

load([dpath,'Example\Ex_results\'...
    '2018-05-18_full_calib_emulator_output_estimates'],...
    'emout');

% Use the true function to find the output at these sample points
y_samps_true = Ex_sim( [2*ones(size(samps_os,1),1) samps_os]);
% Get nondominated outcomes
nondoms = nondominated(y_samps_true);
% Alternative: use GP estimates rather than true function:
nondoms = emout.output_means;
% ( Get direct data from Ex_optimization_workspace, allperfs )
%nondom_ap = nondominated(allperfs);
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');
nondom_ap = ctheta_output_nondom(:,4:6);

% Take a look
Circlesize=50;
%figure; h1 = scatter3(y_samps(:,3),y_samps(:,1),y_samps(:,2),...
%    Circlesize,'b','filled','MarkerFaceAlpha',.4);
figure; h1 = scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    Circlesize,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',.2);
axis vis3d;
hold on;
scatter3(nondoms(:,3),nondoms(:,1),nondoms(:,2),...
    Circlesize,'b','filled','MarkerFaceAlpha',.4);
xlabel('Cost'); ylabel('Deflection'); zlabel('Rotation');
title(['Nondominated direct data (red) and MCMC samples (blue)'...
    ' from full calibration']);

%% Calib parameters post. dist (full calib), compare w/ nondom direct data
clc; clearvars -except dpath ; close all;
% Load MCMC samples
% load([dpath,'Example\Ex_results\'...
%     '2018-05-17_d0_incl_min_cost'],...
%     'results');
load([dpath,'Example\Ex_results\'...
'2018-05-29_set_obs_var_d0'],...
'results');

% Load true samples;
% load([dpath,'Example\Ex_results\'...
%     '2018-05-25_true_ctheta-output'],...
%     'ctheta_output');
% load([dpath,'Example\Ex_results\'...
%     '2018-05-25_true_ctheta-output_nondom'],...
%     'ctheta_output_nondom');
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

figure();
samps = results.samples_os(results.settings.burn_in:end,:);
% Take a look
plot(samps(:,1),samps(:,2),'.');
scatterhist(samps(:,1),samps(:,2),'Kernel','on');
xlim([0,2]); ylim([0,4]);

%%%
% Get allperfs, all_ctheta from Ex_optimization_workspace
%%%

% Find theta settings of nondominated outcomes from direct data
%[nondoms, ndidx] = nondominated(allperfs);
%nondom_theta_true = all_ctheta(ndidx,2:3);
nondom_theta_true = ctheta_output_nondom(:,2:3);

figure();
% Plot theta from nondominated direct outcomes
plot(nondom_theta_true(:,1),nondom_theta_true(:,2),'.');
scatterhist(nondom_theta_true(:,1),nondom_theta_true(:,2),'Kernel','on');
xlim([0,2]); ylim([0,4]); hold on;
scatter1=scatter(samps(:,1),samps(:,2),'.');
scatter1.MarkerFaceAlpha = .1;
scatter1.MarkerEdgeAlpha = .1;

% Get density estimate of direct data and MCMC samples
[ pi_samps , xi_samps , bw_samps ] = ksdensity ( samps );%, 'Bandwidth',...
    %[0.1 0.2]) ; 
[ pi_true  , xi_true  , bw_true ] = ksdensity (nondom_theta_true) ;
% Plot these densities


%% Posterior distribution of calib parameters at specific cost
clc; clearvars -except dpath ; close all;
% Load MCMC samples
load([dpath,'Example\Ex_results\'...
'2018-05-17_cost_grid_results'],...
'results');
% Filter to get all samples with estimated mean near cost_target
cost_target = 20;
dat = [] ; 
for ii = 1:size(results,1)
    newdat = results{ii}.samples_os(results{ii}.settings.burn_in:end,:);
    newdat_cost = results{ii}.model_output.by_sample(:,3);
    newdat = newdat(round(newdat_cost)==cost_target,:);
    dat = [dat ; newdat ] ;
    fprintf('ii: %d, #: %d\n',ii,sum(round(newdat_cost)==cost_target));
end

% take a look
plot(dat(:,1),dat(:,2),'ob');

% Now get the direct data, both unfiltered and nondominated
% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output'],...
    'ctheta_output');
% load([dpath,'Example\Ex_results\'...
%     '2018-05-23_true_ctheta-output_nondom'],...
%     'ctheta_output_nondom');
ctheta_output_at_cost = ...
    ctheta_output(round(ctheta_output(:,6))==cost_target,:);
dd_dat = ctheta_output_at_cost(:,2:3);
% Find nondoms just wrt defl, rotn:
[nd_dd_out, nd_dd_idx] = nondominated(ctheta_output_at_cost(:,4:5)) ; 


nd_dd_dat = ctheta_output_at_cost(nd_dd_idx,2:3);                                                      
%nd_dd_dat = ctheta_output_nondom(...
%    round(ctheta_output_nondom(:,6))==cost_target,2:3);

% Take a look at the direct data (on same plot)
hold on;
plot(dd_dat(:,1),dd_dat(:,2),'og');
plot(nd_dd_dat(:,1),nd_dd_dat(:,2),'or');

%% Plot posterior full calib samps coloring by proximity to pareto front
clc; clearvars -except dpath ; close all;
% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get ranges of outputs
output_ranges = range(ctheta_output(:,4:6));
output_sds    = std(ctheta_output(:,4:6))  ;

% Load MCMC samples (full calib)
% load([dpath,'Example\Ex_results\'...
% '2018-05-28_d0_incl_min_cost'],...
% 'results');
load([dpath,'Example\Ex_results\'...
'2018-05-29_set_obs_var_d0'],...
'results');
samps_os = results.samples_os(results.settings.burn_in:end,:);
% For comparison, also get uniform random samples
unif_samps = ctheta_output(randsample(1:size(ctheta_output,1),...
    size(samps_os,1)),2:3) ; 

% Get true performance at each sample draw
mcmc_perfs = Ex_sim([ 2*ones(size(samps_os,1),1) samps_os ] );
unif_perfs = Ex_sim([ 2*ones(size(unif_samps,1),1) unif_samps ] );

% Get proportions of output ranges
alphas = [0.001 0.01 0.05 0.1 0.5 1 2 3 ];


% Get the true pareto front:
dd_nondoms = ctheta_output_nondom(:,4:6);

% Subtract tol from mcmc_perfs to boost their ``nondominance rating'' by
% tol, and then for each element of mcmc_perfs_boosted, check whether it is
% dominated by any output from the dense direct data grid. (This can be
% done by checking only the nondominated direct data outputs, which are
% stored separately and are loaded above.) If not, then the mcmc sample is
% within tol of the pareto front.
mcmc_perfs_tol_nd_idxs = cell(length(alphas),1) ; % This will store indices
                                                  % of samples within tol
unif_perfs_tol_nd_idxs = cell(length(alphas),1) ; % Similar for unif samps
for jj =1:length(alphas) % Loop through several alpha settings
    alpha = alphas(jj);
    tol = output_sds * alpha ;
%     mcmc_perfs_boosted = mcmc_perfs - tol ;
%     unif_perfs_boosted = unif_perfs - tol ;
    % Get alpha-boosted versions of the outputs
    mcmc_perfs_boosted = [ mcmc_perfs - [tol(1) 0 0 ] ; 
        mcmc_perfs - [0 tol(2) 0 ] ; 
        mcmc_perfs - [0 0 tol(3) ] ] ; 
    unif_perfs_boosted = [ unif_perfs - [tol(1) 0 0 ] ; 
        unif_perfs - [0 tol(2) 0 ] ; 
        unif_perfs - [0 0 tol(3) ] ] ; 
    % Set up for looping through boosted performances
    mcmc_perfs_tol_nd_idx = [] ; % This will be indxs of tol-nondom'd elts
    unif_perfs_tol_nd_idx = [] ;
    for ii = 1:size(mcmc_perfs_boosted,1)
        % samp_ii_nd will be true iff the boosted mcmc samp dominates any
        % element of dd_nondoms (in which case mcmc samp is within tol of
        % the pareto front)
        samp_ii_nd = any(all(mcmc_perfs_boosted(ii,:) < dd_nondoms,2));
        unif_samp_ii_nd = any(all(unif_perfs_boosted(ii,:)<dd_nondoms,2));
        mcmc_perfs_tol_nd_idx = [ mcmc_perfs_tol_nd_idx ;samp_ii_nd ] ;
        unif_perfs_tol_nd_idx = [ unif_perfs_tol_nd_idx ;unif_samp_ii_nd ];
        if mod(ii,1000) == 0
            fprintf('Step %d/%d\n',ii,size(mcmc_perfs_boosted,1));
        end
    end
    % Check outcome: how many samps were within tol
    fprintf('alpha:%d\nMcmc within tol:%d/%d\nUnif within tol:%d/%d\n',...
        alpha,...
        sum(mcmc_perfs_tol_nd_idx),size(mcmc_perfs_boosted,1),...
        sum(unif_perfs_tol_nd_idx),size(unif_perfs_boosted,1));
    mcmc_perfs_tol_nd_idxs{jj} = mcmc_perfs_tol_nd_idx ; 
    unif_perfs_tol_nd_idxs{jj} = unif_perfs_tol_nd_idx ; 
end

samps_within_tol_of_nd.mcmc_perfs_tol_nd_idxs=mcmc_perfs_tol_nd_idxs;
samps_within_tol_of_nd.unif_perfs_tol_nd_idxs=unif_perfs_tol_nd_idxs;
samps_within_tol_of_nd.alphas=alphas;
samps_within_tol_of_nd.output_sds=output_sds;
samps_within_tol_of_nd.mcmc_samps=samps_os;
samps_within_tol_of_nd.unif_samps=unif_samps;
samps_within_tol_of_nd.dd_nondoms=dd_nondoms ;

% load([dpath,'Example\Ex_results\'...
% '2018-05-29_samps_within_tol_of_nd'],...
% 'samps_within_tol_of_nd');

% Find those samples surviving after any of the three perf boosts
mcmc_perfs_tol_nd_any_idxs = cell(size(samps_within_tol_of_nd,1),1);
unif_perfs_tol_nd_any_idxs = cell(size(samps_within_tol_of_nd,1),1);
for ii = 1:size(samps_within_tol_of_nd.alphas,2)
    mx = ...
        samps_within_tol_of_nd.mcmc_perfs_tol_nd_idxs{ii} ; 
    ux = ...
        samps_within_tol_of_nd.unif_perfs_tol_nd_idxs{ii} ; 
    num_samps = size(mx,1)/3;
    mcmc_perfs_tol_nd_any_idx = ...
        any([mx(1:num_samps) mx(num_samps+1:2*num_samps) ...
        mx(2*num_samps+1:end)],2);
    unif_perfs_tol_nd_any_idx = ...
        any([ux(1:num_samps) ux(num_samps+1:2*num_samps) ...
        ux(2*num_samps+1:end)],2);
    fprintf('MCMC within tol:%d/%d\nUnif within tol: %d/%d\n\n',...
        sum(mcmc_perfs_tol_nd_any_idx),num_samps,...
        sum(unif_perfs_tol_nd_any_idx),num_samps) ;
    mcmc_perfs_tol_nd_any_idxs{ii} = mcmc_perfs_tol_nd_any_idx;
    unif_perfs_tol_nd_any_idxs{ii} = unif_perfs_tol_nd_any_idx;
end

% Now make a scatter3 plot with colors indicating alpha level
sz=10 ; % circle size
scatter3(mcmc_perfs(:,1),mcmc_perfs(:,2),mcmc_perfs(:,3),sz); hold on;
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{6},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{6},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{6},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{5},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{5},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{5},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{4},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{4},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{4},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{3},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{3},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{3},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{2},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{2},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{2},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{1},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{1},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{1},3),sz);

%% Plot post. calib samps at spec cost, color by proximity to pareto front
clc; clearvars -except dpath ; close all;
% % Load true samples;
% load([dpath,'Example\Ex_results\'...
%     '2018-05-25_true_ctheta-output'],...
%     'ctheta_output');
% load([dpath,'Example\Ex_results\'...
%     '2018-05-25_true_ctheta-output_nondom'],...
%     'ctheta_output_nondom');
% 
% % Get ranges of outputs
% output_ranges = range(ctheta_output(:,4:6));
% output_sds    = std(ctheta_output(:,4:6))  ;
% 
% Load MCMC samples (full calib)
load([dpath,'Example\Ex_results\'...
'2018-05-17_cost_grid_results'],...
'results');
wh = 5 ; % Which point in grid to use. 5 corresponds to $20.45 des'd cost
samps_os = results{wh}.samples_os(results{wh}.settings.burn_in:end,:);
% % For comparison, also get uniform random samples
% unif_samps = ctheta_output(randsample(1:size(ctheta_output,1),...
%     size(samps_os,1)),2:3) ; 
% 
% Get true performance at each sample draw
mcmc_perfs = Ex_sim([ 2*ones(size(samps_os,1),1) samps_os ] ,true);
% unif_perfs = Ex_sim([ 2*ones(size(unif_samps,1),1) unif_samps ],true );
% 
% % Get proportions of output ranges
% alphas = [0.001 0.01 0.05 0.1 0.5 1 2 3 ];
% 
% 
% % Get the true pareto front:
% dd_nondoms = ctheta_output_nondom(:,4:6);
% 
% % Subtract tol from mcmc_perfs to boost their ``nondominance rating'' by
% % tol, and then for each element of mcmc_perfs_boosted, check whether it is
% % dominated by any output from the dense direct data grid. (This can be
% % done by checking only the nondominated direct data outputs, which are
% % stored separately and are loaded above.) If not, then the mcmc sample is
% % within tol of the pareto front.
% mcmc_perfs_tol_nd_idxs = cell(length(alphas),1) ; % This will store indices
%                                                   % of samples within tol
% unif_perfs_tol_nd_idxs = cell(length(alphas),1) ; % Similar for unif samps
% for jj =1:length(alphas) % Loop through several alpha settings
%     alpha = alphas(jj);
%     tol = output_sds * alpha ;
% %     mcmc_perfs_boosted = mcmc_perfs - tol ;
% %     unif_perfs_boosted = unif_perfs - tol ;
%     % Get alpha-boosted versions of the outputs
%     mcmc_perfs_boosted = [ mcmc_perfs - [tol(1) 0 0 ] ; 
%         mcmc_perfs - [0 tol(2) 0 ] ; 
%         mcmc_perfs - [0 0 tol(3) ] ] ; 
%     unif_perfs_boosted = [ unif_perfs - [tol(1) 0 0 ] ; 
%         unif_perfs - [0 tol(2) 0 ] ; 
%         unif_perfs - [0 0 tol(3) ] ] ; 
%     % Set up for looping through boosted performances
%     mcmc_perfs_tol_nd_idx = [] ; % This will be indxs of tol-nondom'd elts
%     unif_perfs_tol_nd_idx = [] ;
%     for ii = 1:size(mcmc_perfs_boosted,1)
%         % samp_ii_nd will be true iff the boosted mcmc samp dominates any
%         % element of dd_nondoms (in which case mcmc samp is within tol of
%         % the pareto front)
%         samp_ii_nd = any(all(mcmc_perfs_boosted(ii,:) < dd_nondoms,2));
%         unif_samp_ii_nd = any(all(unif_perfs_boosted(ii,:)<dd_nondoms,2));
%         mcmc_perfs_tol_nd_idx = [ mcmc_perfs_tol_nd_idx ;samp_ii_nd ] ;
%         unif_perfs_tol_nd_idx = [ unif_perfs_tol_nd_idx ;unif_samp_ii_nd ];
%         if mod(ii,1000) == 0
%             fprintf('Step %d/%d\n',ii,size(mcmc_perfs_boosted,1));
%         end
%     end
%     % Check outcome: how many samps were within tol
%     fprintf('Mcmc within tol:%d/%d\nUnif within tol:%d/%d\n',...
%         sum(mcmc_perfs_tol_nd_idx),size(mcmc_perfs_boosted,1),...
%         sum(unif_perfs_tol_nd_idx),size(unif_perfs_boosted,1));
%     mcmc_perfs_tol_nd_idxs{jj} = mcmc_perfs_tol_nd_idx ; 
%     unif_perfs_tol_nd_idxs{jj} = unif_perfs_tol_nd_idx ; 
% end
% 
% samps_within_tol_of_nd.mcmc_perfs_tol_nd_idxs=mcmc_perfs_tol_nd_idxs;
% samps_within_tol_of_nd.unif_perfs_tol_nd_idxs=unif_perfs_tol_nd_idxs;
% samps_within_tol_of_nd.alphas=alphas;
% samps_within_tol_of_nd.output_sds=output_sds;
% samps_within_tol_of_nd.mcmc_samps=samps_os;
% samps_within_tol_of_nd.unif_samps=unif_samps;
% samps_within_tol_of_nd.dd_nondoms=dd_nondoms ;
% 
% % save([dpath,'Example\Ex_results\'...
% % '2018-05-29_samps_within_tol_of_nd'],...
% % 'samps_within_tol_of_nd');
% 
% % Find those samples surviving after any of the three perf boosts
% mcmc_perfs_tol_nd_any_idxs = cell(size(samps_within_tol_of_nd,1),1);
% unif_perfs_tol_nd_any_idxs = cell(size(samps_within_tol_of_nd,1),1);
% for ii = 1:size(samps_within_tol_of_nd.alphas,2)
%     mx = ...
%         samps_within_tol_of_nd.mcmc_perfs_tol_nd_idxs{ii} ; 
%     ux = ...
%         samps_within_tol_of_nd.unif_perfs_tol_nd_idxs{ii} ; 
%     num_samps = size(mx,1)/3;
%     mcmc_perfs_tol_nd_any_idx = ...
%         any([mx(1:num_samps) mx(num_samps+1:2*num_samps) ...
%         mx(2*num_samps+1:end)],2);
%     unif_perfs_tol_nd_any_idx = ...
%         any([ux(1:num_samps) ux(num_samps+1:2*num_samps) ...
%         ux(2*num_samps+1:end)],2);
%     fprintf('MCMC within tol:%d/%d\nUnif within tol: %d/%d\n\n',...
%         sum(mcmc_perfs_tol_nd_any_idx),num_samps,...
%         sum(unif_perfs_tol_nd_any_idx),num_samps) ;
%     mcmc_perfs_tol_nd_any_idxs{ii} = mcmc_perfs_tol_nd_any_idx;
%     unif_perfs_tol_nd_any_idxs{ii} = unif_perfs_tol_nd_any_idx;
% end
% 
% samps_within_tol_of_nd.mcmc_perfs_tol_nd_any_idxs = ...
%     mcmc_perfs_tol_nd_any_idxs;
% samps_within_tol_of_nd.unif_perfs_tol_nd_any_idxs = ...
%     unif_perfs_tol_nd_any_idxs;

% All the above was used to create the thing now loaded here:
load([dpath,'Example\Ex_results\'...
'2018-05-30_samps_within_tol_of_nd_with_set_total_obs_var'],...
'samps_within_tol_of_nd');
mcmc_perfs_tol_nd_any_idxs = ...
    samps_within_tol_of_nd.mcmc_perfs_tol_nd_any_idxs;

% Now make a scatter3 plot with colors indicating alpha level
sz=10 ; % circle size
scatter3(mcmc_perfs(:,1),mcmc_perfs(:,2),mcmc_perfs(:,3),sz); hold on;
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{6},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{6},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{6},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{5},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{5},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{5},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{4},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{4},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{4},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{3},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{3},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{3},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{2},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{2},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{2},3),sz);
scatter3(mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{1},1),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{1},2),...
    mcmc_perfs(mcmc_perfs_tol_nd_any_idxs{1},3),sz);

%% Full set total obs var calib with direct data heatmap of proximity to 0
clc; clearvars -except dpath ; close all;
% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get true samples with output closest to 0 (Euclidean distance on
% standardized scale
% First put data on standardized scale
cost_std = (ctheta_output(:,6) - mean(ctheta_output(:,6)))/...
    std(ctheta_output(:,6));
defl_std = (ctheta_output(:,4) - mean(ctheta_output(:,4)))/...
    std(ctheta_output(:,4));
rotn_std = (ctheta_output(:,5) - mean(ctheta_output(:,5)))/...
    std(ctheta_output(:,5));

dd_outputs_std = [defl_std rotn_std cost_std];

% Get zero on standardized scale
zero_pt = -mean(ctheta_output(:,4:6))./std(ctheta_output(:,4:6));

% Now get Euclidean norms of each standardized output
dd_dists = sqrt ( sum ( (dd_outputs_std-zero_pt).^2 , 2 ) ) ;

% Now get the Euclidean norm for all mcmc samples. First, load them:
load([dpath,'Example\Ex_results\'...
'2018-05-29_set_obs_var_d0'],...
'results');
outs = results.model_output.by_sample_true(results.settings.burn_in:end,:);
% Put on standardized scale:
cost_std = (outs(:,3) - mean(ctheta_output(:,6)))/...
    std(ctheta_output(:,6));
defl_std = (outs(:,1) - mean(ctheta_output(:,4)))/...
    std(ctheta_output(:,4));
rotn_std = (outs(:,2) - mean(ctheta_output(:,5)))/...
    std(ctheta_output(:,5));

mcmc_outputs_std = [ defl_std rotn_std cost_std ] ;

% Now get Euclidean norms of each standardized output
mcmc_dists = sqrt ( sum ( (mcmc_outputs_std-zero_pt).^2 , 2 ) ) ;

% Take a look
figure();
scatter(linspace(1,length(mcmc_dists),length(dd_dists)),dd_dists);
hold on;
scatter(1:length(mcmc_dists),mcmc_dists);

% Now take a 3d look at all outputs versus the close direct data outputs
cutoff = quantile(mcmc_dists,.95); % cutoff for close dd output
close_dd_idx = dd_dists <= cutoff; % index of close dd outputs
close_dd_outputs = ctheta_output(close_dd_idx,4:6) ; % close dd outputs
figure();
scatter3(outs(:,1),outs(:,2),outs(:,3),20); axis vis3d; hold on;
scatter3(...
    close_dd_outputs(:,1),close_dd_outputs(:,2),close_dd_outputs(:,3));

% Now take a look at all calib settings at mcmc outputs vs close dd outputs
samps = results.samples_os;
close_dd_theta = ctheta_output(close_dd_idx,2:3);
figure(); scatterhist(samps(:,1),samps(:,2));
figure(); scatterhist(close_dd_theta(:,1),close_dd_theta(:,2));

% Now get a scatterhist of mcmc theta draws with, behind it, all direct
% data theta values colored by Euclidean distance of the standardized
% output to the zero point.
h=figure(); colormap(flipud(jet));
sh=scatterhist(samps(:,1),samps(:,2)); 
hold on; xlim([0 3]); ylim([0 6]);
scatter(ctheta_output(:,2),ctheta_output(:,3),2,dd_dists); hold on;
colorbar('East');
% scatter(samps(:,1),samps(:,2),1,'og','MarkerFaceAlpha',.05,...
%     'MarkerEdgeAlpha',.05);
scatter(samps(:,1),samps(:,2),20,'.g','MarkerFaceAlpha',.5,...
    'MarkerEdgeAlpha',.5);
title({'Posterior \theta draws with marginal distributions'});
xlabel('\theta_1'); ylabel('\theta_2') ;
saveas(h,'FIG_post_theta_settotalvar_w_marginals_and_heatmap.png')

% While we've got that plot up, take a look at the locations of the
% non-dominated thetas.
scatter(ctheta_output_nondom(:,2),ctheta_output_nondom(:,3),'.m');
scatter(samps(:,1),samps(:,2),1,'og','MarkerFaceAlpha',.05,...
    'MarkerEdgeAlpha',.05);
saveas(h,'FIG_post_theta_settotalvar_w_marginals_heatmap_and_nondoms.png');

%% Plot cost grid set total obs var calib with direct data close to 0
clc; clearvars -except dpath ; close all;
% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get true samples w/ defl,rotn output closest to 0 (Euclidean distance on
% standardized scale)
% First put data on standardized scale
cost_std = (ctheta_output(:,6) - mean(ctheta_output(:,6)))/...
    std(ctheta_output(:,6));
defl_std = (ctheta_output(:,4) - mean(ctheta_output(:,4)))/...
    std(ctheta_output(:,4));
rotn_std = (ctheta_output(:,5) - mean(ctheta_output(:,5)))/...
    std(ctheta_output(:,5));

dd_outputs_std = [defl_std rotn_std cost_std];

% Get zero on standardized scale
zero_pt = -mean(ctheta_output(:,4:5))./std(ctheta_output(:,4:5));

% Now get Euclidean norms of each standardized output
dd_dists = sqrt ( sum ( (dd_outputs_std(:,1:2)-zero_pt).^2 , 2 ) ) ;

% Now get the Euclidean norm for all mcmc samples. First, load them:
load([dpath,'Example\Ex_results\'...
'2018-05-30_cost_grid_with_set_total_obs_var'],...
'results');
% Gather all outputs and theta draws together
outs = [];
samps = [];
for ii = 1:size(results,1)
    outs = [outs ; results{ii}.model_output.by_sample_true] ;
    samps = [samps ; results{ii}.samples_os ];
end
% Put on standardized scale:
cost_std = (outs(:,3) - mean(ctheta_output(:,6)))/...
    std(ctheta_output(:,6));
defl_std = (outs(:,1) - mean(ctheta_output(:,4)))/...
    std(ctheta_output(:,4));
rotn_std = (outs(:,2) - mean(ctheta_output(:,5)))/...
    std(ctheta_output(:,5));

mcmc_outputs_std = [ defl_std rotn_std ] ;

% Now get Euclidean norms of each standardized output
mcmc_dists = sqrt ( sum ( (mcmc_outputs_std-zero_pt).^2 , 2 ) ) ;

% Take a look
figure();
scatter(linspace(1,length(mcmc_dists),length(dd_dists)),dd_dists);
hold on;
scatter(1:length(mcmc_dists),mcmc_dists);

% Now take a 3d look at all outputs versus the close direct data outputs
cutoff = quantile(mcmc_dists,.95); % cutoff for close dd output
close_dd_idx = dd_dists <= cutoff; % index of close dd outputs
close_dd_outputs = ctheta_output(close_dd_idx,4:6) ; % close dd outputs
figure();
scatter3(outs(:,1),outs(:,2),outs(:,3),20); axis vis3d; hold on;
scatter3(...
    close_dd_outputs(:,1),close_dd_outputs(:,2),close_dd_outputs(:,3));

% Now take a look at all calib settings at mcmc outputs vs close dd outputs
close_dd_theta = ctheta_output(close_dd_idx,2:3);
figure(); scatterhist(samps(:,1),samps(:,2));
figure(); scatterhist(close_dd_theta(:,1),close_dd_theta(:,2));

% Now get a scatterhist of mcmc theta draws with, behind it, all direct
% data theta values colored by Euclidean distance of the standardized
% output to the zero point.
h=figure(); colormap(flipud(jet));
sh=scatterhist(samps(:,1),samps(:,2),'Kernel','on'); 
hold on; xlim([0 3]); ylim([0 6]);
scatter(ctheta_output(:,2),ctheta_output(:,3),2,dd_dists); hold on;
colorbar('East');
scatter(samps(:,1),samps(:,2),1,'og','MarkerFaceAlpha',.05,...
    'MarkerEdgeAlpha',.05);
title({'Posterior \theta draws with marginal distributions'});
xlabel('\theta_1'); ylabel('\theta_2') ;
%saveas(h,'FIG_post_theta_costgrid_w_marginals_and_heatmap.png')

% While we've got that plot up, take a look at the locations of the
% non-dominated thetas.
scatter(ctheta_output_nondom(:,2),ctheta_output_nondom(:,3),'.m');
scatter(samps(:,1),samps(:,2),1,'og','MarkerFaceAlpha',.05,...
    'MarkerEdgeAlpha',.05);
%saveas(h,'FIG_post_theta_costgrid_w_marginals_heatmap_and_nondoms.png');

% Make a figure like the heatmap above for every point in the cost grid
m = size(results,1);
cost_grid = linspace(15,30,m);
rcosts = round(ctheta_output(:,6)*11/15)*15/11; % Round costs to cost_grid
for ii = 1 : m % Loop through and make figure for each
    ctheta_output_at_cost = ...
        ctheta_output(abs(rcosts - cost_grid(ii))<1e-4,:);
    theta = ctheta_output_at_cost(:,2:3); % Theta vals near cost point
    dd_dists_atcost = dd_dists(abs(rcosts-cost_grid(ii))<1e-4); % Euc dists
    samps = results{ii}.samples_os;
    h = figure(); colormap(flipud(jet));
    sh = scatterhist(samps(:,1),samps(:,2),...
        'Marker','.','Kernel','on');
    hold on; xlim([min(theta(:,1)) max(theta(:,1)) ]);
    ylim([min(theta(:,2)) max(theta(:,2))]);
    scatter(theta(:,1),theta(:,2),500,dd_dists_atcost,'Marker','.'); 
    hold on;
    colorbar('East');
    scatter(samps(:,1),samps(:,2),20,'.g','MarkerFaceAlpha',0.5,...
        'MarkerEdgeAlpha',0.5);
    title(['Posterior \theta draws at cost = '...
        num2str(round(cost_grid(ii),2)) ]);
    xlabel('\theta_1'); ylabel('\theta_2') ;
    saveas(h,['FIG_post_theta_w_marginals_heatmap_atcost_' ...
        num2str(ii) '.png']);
end

% % Make a figure like the heatmap above for every point in the cost grid
% m = size(results,1);
% cost_grid = linspace(15,30,m);
% alpha = 0.25; % This will be the quantiles of mcmc costs used to locate dd
% for ii = 1 : m % Loop through and make figure for each
%     costs = results{ii}.model_output.by_sample_true(:,3); % MCMC costs
%     cq = [quantile(costs,alpha/2) quantile(costs,1-alpha/2)]; % cost quant.
%     cidx =all([ctheta_output(:,6) >= cq(1) ctheta_output(:,6) <= cq(2)],2);
%     ctheta_output_at_cost = ctheta_output(cidx,:);
%     theta = ctheta_output_at_cost(:,2:3); % Theta vals near cost point
%     dd_outs_std_atcost = dd_outputs_std(cidx,:); % Euc dists
%     zero_pt = ([ 0 0 cost_grid(ii) ]-mean(ctheta_output(:,4:6)))...
%         ./std(ctheta_output(:,4:6));
%     dd_dists_atcost = sqrt ( sum ( (dd_outs_std_atcost-zero_pt).^2 , 2 ) );
%     samps = results{ii}.samples_os;
%     h = figure(); colormap(flipud(jet));
%     sh = scatterhist(samps(:,1),samps(:,2),...
%         'Marker','.','Kernel','on');
%     hold on; xlim([0 3 ]);
%     ylim([0 6]);
%     scatter(theta(:,1),theta(:,2),50,...
%         dd_dists_atcost,'Marker','.'); 
%     hold on;
%     colorbar('East');
%     scatter(samps(:,1),samps(:,2),20,'.g','MarkerFaceAlpha',0.5,...
%         'MarkerEdgeAlpha',0.5);
%     title(['Posterior \theta draws at cost = '...
%         num2str(round(cost_grid(ii),2)) ]);
%     xlabel('\theta_1'); ylabel('\theta_2') ;
% end


%% Full prior obs var calib with direct data heatmap of proximity to 0
clc; clearvars -except dpath ; close all;
% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get true samples with output closest to 0 (Euclidean distance on
% standardized scale
% First put data on standardized scale
cost_std = (ctheta_output(:,6) - mean(ctheta_output(:,6)))/...
    std(ctheta_output(:,6));
defl_std = (ctheta_output(:,4) - mean(ctheta_output(:,4)))/...
    std(ctheta_output(:,4));
rotn_std = (ctheta_output(:,5) - mean(ctheta_output(:,5)))/...
    std(ctheta_output(:,5));

dd_outputs_std = [defl_std rotn_std cost_std];

% Get zero on standardized scale
zero_pt = -mean(ctheta_output(:,4:6))./std(ctheta_output(:,4:6));

% Now get Euclidean norms of each standardized output
dd_dists = sqrt ( sum ( (dd_outputs_std-zero_pt).^2 , 2 ) ) ;

% Now get the Euclidean norm for all mcmc samples. First, load them:
load([dpath,'Example\Ex_results\'...
'2018-05-28_d0_incl_min_cost'],...
'results');
outs = results.model_output.by_sample_true(results.settings.burn_in:end,:);
% Put on standardized scale:
cost_std = (outs(:,3) - mean(ctheta_output(:,6)))/...
    std(ctheta_output(:,6));
defl_std = (outs(:,1) - mean(ctheta_output(:,4)))/...
    std(ctheta_output(:,4));
rotn_std = (outs(:,2) - mean(ctheta_output(:,5)))/...
    std(ctheta_output(:,5));

mcmc_outputs_std = [ defl_std rotn_std cost_std ] ;

% Now get Euclidean norms of each standardized output
mcmc_dists = sqrt ( sum ( (mcmc_outputs_std-zero_pt).^2 , 2 ) ) ;

% Take a look
% figure();
% scatter(linspace(1,length(mcmc_dists),length(dd_dists)),dd_dists);
% hold on;
% scatter(1:length(mcmc_dists),mcmc_dists);

% Now take a 3d look at all outputs versus the close direct data outputs
% cutoff = quantile(mcmc_dists,.95); % cutoff for close dd output
% close_dd_idx = dd_dists <= cutoff; % index of close dd outputs
% close_dd_outputs = ctheta_output(close_dd_idx,4:6) ; % close dd outputs
% figure();
% scatter3(outs(:,1),outs(:,2),outs(:,3),20); axis vis3d; hold on;
% scatter3(...
%     close_dd_outputs(:,1),close_dd_outputs(:,2),close_dd_outputs(:,3));
% 
% % Now take a look at all calib settings at mcmc outputs vs close dd outputs
samps = results.samples_os;
% close_dd_theta = ctheta_output(close_dd_idx,2:3);
% figure(); scatterhist(samps(:,1),samps(:,2));
% figure(); scatterhist(close_dd_theta(:,1),close_dd_theta(:,2));

% Now get a scatterhist of mcmc theta draws with, behind it, all direct
% data theta values colored by Euclidean distance of the standardized
% output to the zero point.
h=figure(); colormap(flipud(jet));
sh=scatterhist(samps(:,1),samps(:,2)); 
hold on; xlim([0 3]); ylim([0 6]);
scatter(ctheta_output(:,2),ctheta_output(:,3),2,dd_dists); hold on;
colorbar('East');
scatter(samps(:,1),samps(:,2),20,'.g','MarkerFaceAlpha',.5,...
    'MarkerEdgeAlpha',.5);
title({'Posterior \theta draws with marginal distributions'});
xlabel('\theta_1'); ylabel('\theta_2') ;
saveas(h,'FIG_hmfc.png')

% While we've got that plot up, take a look at the locations of the
% non-dominated thetas.
scatter(ctheta_output_nondom(:,2),ctheta_output_nondom(:,3),'.m');
scatter(samps(:,1),samps(:,2),20,'.g','MarkerFaceAlpha',.5,...
    'MarkerEdgeAlpha',.5);
saveas(h,'FIG_hmfc_nd.png');

%% Show mahalanobis distance of samples from nondominated pts
% Both from the calibration using set total observation variance, and the
% older version using a prior on observation variance.
clc; clearvars -except dpath ; close all;
% Load true samples;
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output'],...
    'ctheta_output');
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

load([dpath,'Example\Ex_results\'...
    '2018-05-28_d0_incl_min_cost'],...
    'results');

% Get nondominated theta points by themselves
nondoms = ctheta_output_nondom(:,2:3);

% Get sample draws
samps = results.samples_os(results.settings.burn_in+2:end,:);

% Get random (uniform) sample from the parameter space
unif_samps = rand(size(samps)) .* [3 6];

% Now get info o mahalanobis distance of unif_samps and samps from nondoms
unif_md = mahal(unif_samps,nondoms);
samps_md = mahal(samps,nondoms);
n=size(unif_md,1);
%plot(1:n,unif_md,'.'); hold on; plot(1:n,samps_md,'.');

% Get histograms of mahalanobis distance of samples from nondoms
f=figure();
h=histogram(samps_md);
title('Mahalanobis distance of samples from nondominated region');
xlabel('Distance'); ylabel('Number of samples');
ylim([0 4000]);
saveas(h,'FIG_md.png');
hold on;
title({'Mahalanobis distance of samples from nondominated region:'...
    'comparison with uniformly sampled points (red)'});
histogram(unif_md);
saveas(h,'FIG_mdwus.png');
% 
% % Get scatterhist of samps with nondoms
% figure();
% h=scatterhist(samps(:,1),samps(:,2),'Marker','.'); hold on;
% scatter(nondoms(:,1),nondoms(:,2),'Marker','.','MarkerFaceAlpha',0.5);
% scatter(samps(:,1),samps(:,2),15,[0 0.447 .741],'Marker','.',...
%     'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.4);
% ttl=title('Posterior distribution with marginals and nondominated region');
% xlabel('\theta_1');ylabel('\theta_2');
% ttl.Position=[.8 6.5 0];
% saveas(h,'FIG_post_dist_w_marginals_and_nondoms_for_prior_obs_var.png');
% 
% % Now do the above stuff for the new calib version, set total obs var
% load([dpath,'Example\Ex_results\'...
%     '2018-05-29_set_obs_var_d0'],...
%     'results');
% samps = results.samples_os(results.settings.burn_in+2:end,:);
% samps_md = mahal(samps,nondoms);
% figure();
% h=scatterhist(samps(:,1),samps(:,2),'Marker','.'); hold on;
% xlim([0,3]);ylim([0,6]);
% scatter(nondoms(:,1),nondoms(:,2),'Marker','.','MarkerFaceAlpha',0.5);
% scatter(samps(:,1),samps(:,2),15,[0 0.447 .741],'Marker','.',...
%     'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.4);
% ttl=title({['Posterior distribution with marginals '...
%     'and nondominated region,'],'using set total observation variance' });
% xlabel('\theta_1');ylabel('\theta_2');
% ttl.Position=[.8 6.2 0];
% saveas(h,'FIG_post_dist_w_marginals_and_nondoms_for_prior_obs_var.png');

% Now get the same hist of mahalanobis distances for the new method samps
load([dpath,'Example\Ex_results\'...
    '2018-05-29_set_obs_var_d0'],...
    'results');
samps = results.samples_os(results.settings.burn_in+2:end,:);
samps_md = mahal(samps,nondoms);
% Get histograms of mahalanobis distance of samples from nondoms
figure();
h=histogram(samps_md);
h.BinWidth=0.35; ylim([0 4000]);
title('Mahalanobis distance of samples from nondominated region');
xlabel('Distance'); ylabel('Number of samples');
saveas(h,'FIG_md_stov.png');
hold on;
histogram(unif_md);
title({'Mahalanobis distance of samples from nondominated region'...
    'comparison with uniformly sampled points (red)'});
saveas(h,'FIG_mdwus_stov.png');

% Now compare samps to the nondoms themselves
figure(); h=histogram(samps_md);
h.BinWidth=0.35; ylim([0 4000]);
title({'Mahalanobis distance of samples from nondominated region:'...
    'comparison to random sample of nondominated points (red)'});
xlabel('Distance'); ylabel('Number of samples');
hold on;
nondom_samp = nondoms(randsample(1:size(nondoms,1),size(samps,1)),:);
nd_md = mahal(nondom_samp,nondom_samp);
histogram(nd_md,'BinWidth',0.35);
saveas(h,'FIG_mdwnd_stov.png');
