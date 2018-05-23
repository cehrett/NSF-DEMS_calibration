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

%% Get deflection, rotation "bands" as function of cost, with true trade-off
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
data_true = nondominated(allperfs); % Get allperfs 
                               % from Ex_optimization_workspace
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
load([dpath,'Example\Ex_results\'...
'2018-05-17_d0_incl_min_cost'],...
'results'); % This is the results
load([dpath,'Example\Ex_results\'...
    '2018-05-18_full_calib_emulator_output_estimates'],...
    'emout'); % This is the GP estimate of output at each draw of the MCMC

cost_target = 20;

% Get outputs at that cost
samps_at_cost_idx = round(emout(:,3) == cost_target);
samps_at_cost = emout(samps_at_cost_idx,:);

% Take a look
plot(samps_at_cost(:,1),samps_at_cost(:,2),'o');

% Now find the same using direct data
load([dpath,'Example\Ex_results\'...
    '2018-05-23_true_ctheta-output'],...
    'ctheta_output');
dd_at_cost_idx = round(ctheta_output(:,6) == cost_target);
dd_at_cost = ctheta_output(:,4:5);

% Take a look at direct data tradeoff
hold on;
plot(dd_at_cost(:,1),dd_at_cost(:,2),'o');


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
nondoms = nondominated(emout.output_means);
% ( Get direct data from Ex_optimization_workspace, allperfs )
%nondom_ap = nondominated(allperfs);
load([dpath,'Example\Ex_results\'...
    '2018-05-18_direct_data_nondom'],...
    'nondom_ap');

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

% Compare with data obtained directly
scatter3(nondom_ap(:,3),nondom_ap(:,1),nondom_ap(:,2),...
    Circlesize,nondom_ap(:,3),'r','filled','MarkerFaceAlpha',.2);

%% Get 3d figure of nondominated direct data plus nondom cost grid MCMC
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
    '2018-05-18_direct_data_nondom'],...
    'nondom_ap');

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
    '2018-05-18_direct_data_nondom'],...
    'nondom_ap');

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
load([dpath,'Example\Ex_results\'...
    '2018-05-17_d0_incl_min_cost'],...
    'results');

% Load true samples;
ctheta_output = [all_ctheta allperfs];
save([dpath,'Example\Ex_results\'...
    '2018-05-23_true_ctheta-output'],...
    'ctheta_output');

samps = results.samples_os(results.settings.burn_in:end,:);
% Take a look
plot(samps(:,1),samps(:,2),'.');
scatterhist(samps(:,1),samps(:,2),'Kernel','on');

%%%
% Get allperfs, all_ctheta from Ex_optimization_workspace
%%%

% Find theta settings of nondominated outcomes from direct data
[nondoms, ndidx] = nondominated(allperfs);
nondom_theta_true = all_ctheta(ndidx,2:3);

% Plot theta from nondominated direct outcomes
plot(nondom_theta_true(:,1),nondom_theta_true(:,2),'.');
scatterhist(nondom_theta_true(:,1),nondom_theta_true(:,2),'Kernel','on');

% Get density estimate of direct data and MCMC samples
[ pi_samps , xi_samps ] = ksdensity ( samps ) ; 
[ pi_true  , xi_true  ] = ksdensity (nondom_theta_true) ;
% Plot these densities


