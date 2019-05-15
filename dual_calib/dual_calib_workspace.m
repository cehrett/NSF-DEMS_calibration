% Dual calibration workspace
%% Set path string and add paths
clc; clear all; close all;

direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'dual_calib']);
addpath([dpath,'Example']);

%% Explore output of example function
clc ; clearvars -except dpath ; close all ;

% Define inputs
xmin = .5;
xrange = .5;
x = linspace(0,1);
t1min = 1.5;
t1range = 3;
t1=linspace(0,1);
t2min = 0;
t2range = 5;
t2 = linspace(0,1);

[X,T1,T2] = meshgrid(x,t1,t2) ; 
Y = reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1),length(x),length(t1),length(t2));

% Take a look
xidx=100;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
surf(tt1 * t1range + t1min,tt2*t2range+t2min,reshape(Y(:,xidx,:),100,100));

%% Gather data to use for emulator and for "real" observations
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Get simulation design points and observations
X = lhsdesign(250,3);
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;
% Get simulation observations
sim_y = dual_calib_example_fn(sim_x,sim_t1,sim_t2);

% Now set "real" theta1 and get "real" observations
theta1 = 2;
X = lhsdesign(30,2);
obs_x = X(:,1) * xrange + xmin;
obs_t2 = X(:,2) * t2range + t2min;
% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;
% Get "real" observations without noise or discrepancy
obs_y_noiseless = dual_calib_example_fn(obs_x,obs_t1,obs_t2);

% Now noise it up (still keep no discrepancy for now)
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
std_y = std(sim_y);
mean_y = mean(sim_y);
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;

% Now set desired observations
des_x = linspace(0,1,8)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Now package everything up and save it
clear X sigma obs_y_noiseless dpath
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-01-15_dual_calib_raw_data']);
% And get dpath back:
direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 


%% Get MLEs for covariance parameters for emulator GP
clc ; clearvars -except dpath ; close all ; 

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-01-15_dual_calib_raw_data']);

% Normalize the simulation inputs, standardize the outputs
x = (sim_x - xmin) ./ xrange ;
t1 = (sim_t1 - t1min) ./ t1range;
t2 = (sim_t2 - t2min) ./ t2range;
y = (sim_y - mean_y) ./ std_y;

% Define function for minimization
f = @(rl) ...
    -logmvnpdf(y',0,...
    gp_cov([rl(1) rl(2) rl(3)],[x t1 t2],[x t1 t2],rl(4),false) + ...
    1e-4*eye(size(x,1)));

% Perform minimization
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0 0 0];
ub = [1 1 1 Inf];
x0 = [.5 .5 .5 1];
[inp,fval,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,lb,ub);

% Result: rho values 0.9929, 0.7844, 0.0779; lambda value 0.0836.


%% Perform dual calibration without emulator or discrepancy
% Here calibration will be performed without emulator or (true) discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-01-15_dual_calib_raw_data']);
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'EmulatorCovHypers',[.5 .5 .5 Inf],'obs_discrep',false);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-01-25_dual_calib_nosim_nodiscr'],'results');

%% Perform dual calibration with emulator and without discrepancy
% Here calibration will be performed with emulator but without discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-01-15_dual_calib_raw_data']);

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'obs_discrep',false);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-08_dual_calib_nodiscr'],'results');

%% Gather data for emulator and for "real" observations, with discrepancy
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Get simulation design points and observations
X = lhsdesign(250,3);
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;
% Get simulation observations
sim_y = dual_calib_example_fn(sim_x,sim_t1,sim_t2);

% Now set "real" theta1 and get "real" observations
theta1 = 2;
X = lhsdesign(30,2);
obs_x = X(:,1) * xrange + xmin;
obs_t2 = X(:,2) * t2range + t2min;
% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;
% Get "real" observations without noise but with discrepancy
discrep=2;
obs_y_noiseless = dual_calib_example_fn(X(:,1),xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,X(:,2),t2min,t2range,0,1,discrep);

% Now noise it up (still keep no discrepancy for now)
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
std_y = std(sim_y);
mean_y = mean(sim_y);
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;

% Now set desired observations
des_x = linspace(0,1,8)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Now package everything up and save it
clear X sigma obs_y_noiseless dpath
% load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-14_dual_calib_raw_data']);
% And get dpath back:
direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

%% Perform dual calibration without emulator but with discrepancy
% Here calibration will be performed without emulator or (true) discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-02-14_dual_calib_raw_data']);
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'EmulatorCovHypers',[.5 .5 .5 Inf],'obs_discrep',true,...
    'emulator',false,'modular',false);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-14_dual_calib_nosim_discr'],'results');

%% Perform dual calibration with emulator and discrepancy
% Here calibration will be performed with emulator but without discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-02-14_dual_calib_raw_data']);

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'obs_discrep',true);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-14_dual_calib'],'results');

%% Get data with discrepancy, perform dual calibration, and plot results
clc ; clearvars -except dpath ; close all ;

% Set real theta1, discrepancy version, whether modular
theta1 = 2;
discrep = 1;
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Commenting out the following lines because design is now saved
% % Get simulation design points and observations
% X = lhsdesign(250,3);
% sim_x = X(:,1) * xrange + xmin;
% sim_t1 = X(:,2) * t1range + t1min;
% sim_t2 = X(:,3) * t2range + t2min;
% % Get simulation observations
% sim_y = dual_calib_example_fn(sim_x,sim_t1,sim_t2);
% 
% % Now get "real" observations
% X = lhsdesign(50,2);
% obs_x = X(:,1) * xrange + xmin;
% obs_t2 = X(:,2) * t2range + t2min;
% std_y = std(sim_y);
% mean_y = mean(sim_y);

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2; mean_y = design.mean_y;
std_y = design.std_y;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;
% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;

% Now set desired observations
des_x = linspace(0,1,8)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t)...
    dual_calib_example_fn(x,0,1,theta1*ones(size(x))...
    ,0,1,t,0,1,0,1,discrep)-...
    dual_calib_example_fn(x,0,1,theta1*ones(size(x)),0,1,t,0,1,...
    mean_y,std_y,0);
avg_disc = integral2(int_fn,xmin,xmin+xrange,t2min,t2min+t2range) / ...
    (xrange * t2range) ; 
avg_disc=0 ; 
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 
% The below lines are included so that an informative discrepancy mean can
% be set for discrep=3. 
% c = 0.055 ; % This controls the size of the additive discrepancy
% d = 0; %Raising/lowering c makes discrepancy more/less aggressive
% obs_discrep_mean = @(xx,tt2) ...
%     (c*(xx*xrange+xmin).*(tt2*t2range+t2min)+d)/std_y;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',1e4,'burn_in',.4,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
results = MCMC_dual_calib(settings);

burn_in = results.settings.burn_in;
min_t1 = results.settings.min_t1; range_t1 = results.settings.range_t1;
min_t2 = results.settings.min_t2; range_t2 = results.settings.range_t2;

% First, get prior and posterior theta1
f1 = figure('pos',[10 10 400 300]);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','Posterior dist.','True value');
title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
% Save it:
set(f1,'color','white');
% saveas(f1,'FIG_dual_calib_post_theta1-2.png');

% Second, get prior and posterior theta2
f2 = figure('pos',[420 10 400 300]);
left_color = [.5 .5 0];
right_color = [.5 0 .5];
set(f2,'defaultAxesColorOrder',[left_color; right_color]);
% Plot prior
fill([min_t2 min_t2 + range_t2 min_t2 + range_t2 min_t2],...
    [0 0 1/range_t2 1/range_t2],'g','EdgeColor','none');
xlim([min_t2 min_t2 + range_t2]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85);
% Get and plot true theta2
fmfn = @(z) dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],0,5);
yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','Posterior dist.','Optimal value');
title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');

% Also examine discrepancy
% Define inputs
xmin = .5;
xrange = .5;
x = linspace(0,1);
t1min = 1.5;
t1range = 3;
t1=linspace(0,1);
t2min = 0;
t2range = 5;
t2 = linspace(0,1);

[X,T1,T2] = meshgrid(x,t1,t2) ; 
Y = reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
Yd= reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

% Take a look
f3=figure('pos',[680 140 600 500]);
xidx=50;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
Discrep = Yd-Y;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Y(:,xidx,:),100,100),'EdgeAlpha',.25);
hold on;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Yd(:,xidx,:),100,100),'EdgeAlpha',.25);
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Discrep(:,xidx,:),100,100),'EdgeAlpha',.25);

% Save results
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',int2str(discrep)]);
% save(locstr,'results');


%% Perform KOH calibration, and plot results
clc ; clearvars -except dpath ; close all ;

% Set real theta1, discrepancy version, whether modular
theta1 = 2;
discrep = 0;
modular = false;
obs_discrep = false; % Whether to include discrep term for real obs

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2; mean_y = design.mean_y;
std_y = design.std_y;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations
% Get same observations used in DCTO version:
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',int2str(discrep)]);
load(locstr,'results');
obs_y = results.settings.obs_y * std_y + mean_y;

% Now set desired observations; since we're doing KOH, there are none
des_x = []; des_y = [];

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t)...
    dual_calib_example_fn(x,0,1,theta1*ones(size(x))...
    ,0,1,t,0,1,0,1,discrep)-...
    dual_calib_example_fn(x,0,1,theta1*ones(size(x)),0,1,t,0,1,...
    mean_y,std_y,0);
avg_disc = integral2(int_fn,xmin,xmin+xrange,t2min,t2min+t2range) / ...
    (xrange * t2range) ; 
avg_disc=0 ; 
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 
% The below lines are included so that an informative discrepancy mean can
% be set for discrep=3. 
% c = 0.055 ; % This controls the size of the additive discrepancy
% d = 0; %Raising/lowering c makes discrepancy more/less aggressive
% obs_discrep_mean = @(xx,tt2) ...
%     (c*(xx*xrange+xmin).*(tt2*t2range+t2min)+d)/std_y;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',1e4,'burn_in',.4,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',modular);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
results = MCMC_dual_calib(settings);

burn_in = results.settings.burn_in;
min_t1 = results.settings.min_t1; range_t1 = results.settings.range_t1;
min_t2 = results.settings.min_t2; range_t2 = results.settings.range_t2;

% First, get prior and posterior theta1
f1 = figure('pos',[10 10 400 300]);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','Posterior dist.','True value');
title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
% Save it:
set(f1,'color','white');
% saveas(f1,'FIG_dual_calib_post_theta1-2.png');

% Second, get prior and posterior theta2
f2 = figure('pos',[420 10 400 300]);
left_color = [.5 .5 0];
right_color = [.5 0 .5];
set(f2,'defaultAxesColorOrder',[left_color; right_color]);
% Plot prior
fill([min_t2 min_t2 + range_t2 min_t2 + range_t2 min_t2],...
    [0 0 1/range_t2 1/range_t2],'g','EdgeColor','none');
xlim([min_t2 min_t2 + range_t2]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85);
% Get and plot true theta2
fmfn = @(z) dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],0,5);
yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','Posterior dist.','Optimal value');
title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');

% Also examine discrepancy
% Define inputs
xmin = .5;
xrange = .5;
x = linspace(0,1);
t1min = 1.5;
t1range = 3;
t1=linspace(0,1);
t2min = 0;
t2range = 5;
t2 = linspace(0,1);

[X,T1,T2] = meshgrid(x,t1,t2) ; 
Y = reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
Yd= reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

% Take a look
f3=figure('pos',[830 10 400 300]);
xidx=50;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
Discrep = Yd-Y;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Y(:,xidx,:),100,100),'EdgeAlpha',.25);
hold on;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Yd(:,xidx,:),100,100),'EdgeAlpha',.25);
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Discrep(:,xidx,:),100,100),'EdgeAlpha',.25);

% Save results
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_KOH_calib_discrep',int2str(discrep)]);
% save(locstr,'results');

%% Perform CTO after KOH calibration, compare to DCTO
clc ; clearvars -except dpath ; close all ; 

% Which discrepancy version are we working with:
discrep = 5;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;
theta1=2;

% Load saved response mean and std
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
mean_y = design.mean_y; std_y = design.std_y;
clear design ;

% Load results of KOH calibration of theta1
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_KOH_calib_discrep',int2str(discrep)]);
load(locstr,'results');

theta1_samps = ...
    (results.theta1(results.settings.burn_in:end) - ...
    results.settings.min_t1) / results.settings.range_t1 ;
clear results ; % Don't need these any more

% Define objective function
mean_sim = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,randsample(theta1_samps,1),t1min,t1range,...
    b,t2min,t2range,...
    mean_y,std_y);

% Now set desired observations
obs_x = linspace(0,1,8)' * xrange + xmin;
obs_y = zeros(size(obs_x,1),1);
obs_t2 = obs_x ; % obs_t2 doesn't matter, it is ignored by the obj fn
% Notice we named the above obs_x,obs_y rather than des_x,des_y. This is
% because we are performing KOH calibration on theta2. We set des_x,des_y
% to be empty, since we only have a KOH calibration here, no DCTO secondary
% calibration.
des_x = [] ; des_y = [] ;

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
avg_disc=0 ; 
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 


% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',2e4,'burn_in',.2,...
    'obs_discrep',true,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
results = MCMC_dual_calib(settings);

%%% Now make histograms of results, including DCTO results
burn_in = results.settings.burn_in;
min_t1 = t1min; range_t1 = t1range;
min_t2 = t2min; range_t2 = t2range;
theta1_samps = theta1_samps * range_t1 + min_t1; % Return t1 to orig scale

% Load DCTO results
KOH_results = results;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',int2str(discrep)]);
load(locstr,'results');

% First, get prior and posterior theta1
f1 = figure('pos',[10 10 400 300]);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(theta1_samps,'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
histogram(results.theta1(burn_in:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
% Save it:
set(f1,'color','white');
% saveas(f1,'FIG_dual_calib_post_theta1-2.png');

% Second, get prior and posterior theta2
f2 = figure('pos',[420 10 400 300]);
left_color = [.5 .5 0];
right_color = [.5 0 .5];
set(f2,'defaultAxesColorOrder',[left_color; right_color]);
% Plot prior
fill([min_t2 min_t2 + range_t2 min_t2 + range_t2 min_t2],...
    [0 0 1/range_t2 1/range_t2],'g','EdgeColor','none');
xlim([min_t2 min_t2 + range_t2]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(KOH_results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
fmfn = @(z) dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','KOH','DCTO','Optimal value');
title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');


% Save results
% results = KOH_results;
% locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-03-29_KOH_calib_step2_discrep',int2str(discrep),'_inf']);
% save(locstr,'results');

%% Find MSEs of DCTO and KOH for each of 7 versions
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Set true theta1
theta1 = 2;

% Set up containers for MSEs and vars
t1_mse = zeros(2,8);
t2_mse = zeros(2,8);
t1_var = zeros(2,8);
t2_var = zeros(2,8);

% Loop through each discrep value and find MSEs
for ii=0:7
    discrep = ii;
    
    % Set up to load appropriate files
    discrep_str = int2str(discrep);
    if ii==7 discrep = 5; discrep_str = [int2str(discrep),'_inf']; end
    
    % Find optimal theta2
    fmfn = @(z)dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
    theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
    
    % Load DCTO and KOH results
    dctostr = sprintf(['C:\\Users\\carle\\Documents',...
        '\\MATLAB\\NSF DEMS\\Phase 1\\',...
        'dual_calib\\dual_calib_stored_data\\'...
        '2019-03-28_DCTO_discrep',discrep_str]);
    koh1str = sprintf(['C:\\Users\\carle\\Documents',...
        '\\MATLAB\\NSF DEMS\\Phase 1\\',...
        'dual_calib\\dual_calib_stored_data\\'...
        '2019-03-28_KOH_calib_discrep',discrep_str]);
    koh2str = sprintf(['C:\\Users\\carle\\Documents',...
        '\\MATLAB\\NSF DEMS\\Phase 1\\',...
        'dual_calib\\dual_calib_stored_data\\'...
        '2019-03-29_KOH_calib_step2_discrep',discrep_str]);
    load(dctostr); 
    burnin=results.settings.burn_in; 
    % Get theta1 and theta2 from dcto
    dcto_t1 = results.theta1(burnin:end); n = length(dcto_t1);
    dcto_t2 = results.theta2(burnin:end);
    % Get theta1 and theta2 from koh and cto
    load(koh1str); koh_t1 = results.theta1(burnin:end);
    load(koh2str); koh_t2 = results.theta1(burnin:end);
    
    % Find MSEs
    dcto_t1_mse = sum((dcto_t1 - theta1).^2)/n;
    dcto_t2_mse = sum((dcto_t2 - theta2).^2)/n;
    koh_t1_mse = sum((koh_t1 - theta1).^2)/n;
    koh_t2_mse = sum((koh_t2 - theta2).^2)/n;
    t1_mse(:,ii+1) = [dcto_t1_mse koh_t1_mse];
    t2_mse(:,ii+1) = [dcto_t2_mse koh_t2_mse];
    
    % Also var
    dcto_t1_var = var(dcto_t1) ; 
    dcto_t2_var = var(dcto_t2) ;    
    koh_t1_var = var(koh_t1) ;
    koh_t2_var = var(koh_t2) ;
    t1_var(:,ii+1) = [dcto_t1_var koh_t1_var] ;
    t2_var(:,ii+1) = [dcto_t2_var koh_t2_var] ;
    
end

t1_mse
t2_mse
t1_var
t2_var

%% Take a look at realizations of discrepancies
clc ; clearvars -except dpath ; close all ;

%%% Open figure and set which rho,lambda draws are used
f=figure('pos',[680 140 600 500]);
idx = randsample(4000,1);

%%% Set which discrepancy version is being examined
discrep = 7; 

%%% Define input rescaling settings
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;
theta1 = 2;

%%% Load the MCMC results
discrep_str = int2str(discrep);
if discrep == 7 discrep_str = ['5_inf'] ; discrep = 5 ; end
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',discrep_str]);
load(locstr);
burnin = results.settings.burn_in;
rhos = results.obs_rho(burnin:end,:) ;
lambdas = results.obs_lambda(burnin:end,:);
theta1s = results.theta1(burnin:end,:);
n = size(rhos,1);

%%% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-4*eye(size(X)); % adds nugget for computational stablty
prior_mean = results.settings.mean_obs;
prior_cov = @(rho,x,t2,xp,t2p,lambda) ... 
    gp_cov(rho,[x t2],[xp t2p],lambda,false);
updated_mean = @(y,x,t2,xnew,t2new,rho,lambda) ...
    prior_mean(xnew,t2new) + ...
    add_nug(prior_cov(rho,xnew,t2new,x,t2,lambda)) * ...
    inv(add_nug(prior_cov(rho,x,t2,x,t2,lambda))) * ...
    (y - prior_mean(x,t2)) ;
updated_cov = @(x,t2,xnew,t2new,rho,lambda) ...
    add_nug(prior_cov(rho,xnew,t2new,xnew,t2new,lambda)) - ...
    add_nug(prior_cov(rho,xnew,t2new,x,t2,lambda)) * ...
    inv(add_nug(prior_cov(rho,x,t2,x,t2,lambda))) * ...
    add_nug(prior_cov(rho,x,t2,xnew,t2new,lambda)) ;

%%% Get vectors for observed points and prediction points
obs_x = results.settings.obs_x ; 
obs_t2 = results.settings.obs_t2 ; 
obs_y = results.settings.obs_y ; 
h=40; % Grid mesh setting for theta2
[xp, t2p] = meshgrid(linspace(0,1,5),linspace(0,1,h)) ;
xp=xp(:) ; t2p = t2p(:);

%%% Get mean of discrep for a random draw of rho,lambda
rho = rhos(idx,:) ; lambda = lambdas(idx,:) ; theta1 = theta1s(idx,:) ;
obs_disc = obs_y - dual_calib_example_fn(obs_x,xmin,xrange,...
    theta1,0,1,obs_t2,t2min,t2range,...
    results.settings.mean_y,results.settings.std_y,0);
% Gather discrepancy mean on standardized scale:
d_std = updated_mean(obs_disc,obs_x,obs_t2,xp,t2p,rho,lambda) ;
% Transform to original scale:
d = d_std * results.settings.std_y ;
% Now get standard deviations of d:
d_std_cov = updated_cov(obs_x,obs_t2,xp,t2p,rho,lambda) ; 
d_std_sds = sqrt(diag(d_std_cov)+0.05) ; 
d_sds = d_std_sds * results.settings.std_y;

%%% Plot the discrepancy
% Define inputs
x = 0.5;
t1=linspace(0,1);
t2 = linspace(0,1);

[X,T1,T2] = meshgrid(x,t1,t2) ; 
Y = reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
Yd= reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

% Take a look
xx=reshape(X,100,100);
tt1=reshape(T1,100,100);
tt2=reshape(T2,100,100);
Discrep = Yd-Y;
% surf(tt1*t1range+t1min,tt2*t2range+t2min,...
%     reshape(Y,100,100),'EdgeAlpha',.25);
% hold on;
% surf(tt1*t1range+t1min,tt2*t2range+t2min,...
%     reshape(Yd,100,100),'EdgeAlpha',.25);
% True discrepancy:
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Discrep,100,100),'EdgeAlpha',.25);
hold on;
% Now we add the estimated discrepancy:
[X,T1,T2] = meshgrid(x,t1,t2p(1:size(unique(t2p),1))) ; 
xx=reshape(X,100,h);
tt1=reshape(T1,100,h);
tt2=reshape(T2,100,h);
dd = repmat(d(xp==0.5)',100,1);
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    dd,'EdgeAlpha',.25);
% Now add upper and lower 2sd surfaces
d_upper = d + 2 * d_sds ; d_lower = d - 2 * d_sds ;
dd_upper = repmat(d_upper(xp==.5)',100,1);
dd_lower = repmat(d_lower(xp==.5)',100,1);
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    dd_upper,'EdgeAlpha',.25);
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    dd_lower,'EdgeAlpha',.25);

f.Children.View = [-90 0]; hold off ; pause(0.01);

%% Perform DCTO with informatively chosen target outcomes
clc ; clearvars -except dpath ; close all ;

% Set real theta1, optimal theta2, discrepancy version, whether modular
theta1 = 2;
discrep = 6;
if discrep < 5 theta2=4/3 ; else theta2=1 ; end
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2; mean_y = design.mean_y;
std_y = design.std_y;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations
% Get same observations used in DCTO version:
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',int2str(discrep)]);
load(locstr,'results');
obs_y = results.settings.obs_y * std_y + mean_y;

% Now set desired observations, informatively
des_x = linspace(0,1,8)' * xrange + xmin;
des_y = dual_calib_example_fn(des_x,theta1,theta2) - 0.01;%des_y=des_y*0;

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t)...
    dual_calib_example_fn(x,0,1,theta1*ones(size(x))...
    ,0,1,t,0,1,0,1,discrep)-...
    dual_calib_example_fn(x,0,1,theta1*ones(size(x)),0,1,t,0,1,...
    mean_y,std_y,0);
avg_disc = integral2(int_fn,xmin,xmin+xrange,t2min,t2min+t2range) / ...
    (xrange * t2range) ; 
avg_disc=0 ; 
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 
% The below lines are included so that an informative discrepancy mean can
% be set for discrep=3. 
% c = 0.055 ; % This controls the size of the additive discrepancy
% d = 0; %Raising/lowering c makes discrepancy more/less aggressive
% obs_discrep_mean = @(xx,tt2) ...
%     (c*(xx*xrange+xmin).*(tt2*t2range+t2min)+d)/std_y;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',1e4,'burn_in',.4,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular);

% Perform calibration
results = MCMC_dual_calib(settings);

burn_in = results.settings.burn_in;
min_t1 = results.settings.min_t1; range_t1 = results.settings.range_t1;
min_t2 = results.settings.min_t2; range_t2 = results.settings.range_t2;

% First, get prior and posterior theta1
f1 = figure('pos',[10 10 400 300]);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','Posterior dist.','True value');
title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
% Save it:
set(f1,'color','white');
% saveas(f1,'FIG_dual_calib_post_theta1-2.png');

% Second, get prior and posterior theta2
f2 = figure('pos',[420 10 400 300]);
left_color = [.5 .5 0];
right_color = [.5 0 .5];
set(f2,'defaultAxesColorOrder',[left_color; right_color]);
% Plot prior
fill([min_t2 min_t2 + range_t2 min_t2 + range_t2 min_t2],...
    [0 0 1/range_t2 1/range_t2],'g','EdgeColor','none');
xlim([min_t2 min_t2 + range_t2]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85);
% Get and plot true theta2
fmfn = @(z) dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],0,5);
yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','Posterior dist.','Optimal value');
title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');

% Also examine discrepancy
% Define inputs
xmin = .5;
xrange = .5;
x = linspace(0,1);
t1min = 1.5;
t1range = 3;
t1=linspace(0,1);
t2min = 0;
t2range = 5;
t2 = linspace(0,1);

[X,T1,T2] = meshgrid(x,t1,t2) ; 
Y = reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
Yd= reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

% Take a look
f3=figure('pos',[680 140 600 500]);
xidx=50;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
Discrep = Yd-Y;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Y(:,xidx,:),100,100),'EdgeAlpha',.25);
hold on;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Yd(:,xidx,:),100,100),'EdgeAlpha',.25);
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Discrep(:,xidx,:),100,100),'EdgeAlpha',.25);

% Check posterior means
disp(mean(results.theta1));
disp(mean(results.theta2));

% Save results
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-05-01_DCTO_discrep',int2str(discrep)]);
% save(locstr,'results');

%% Get predictive distribution and compare to Pareto front
clc ; clearvars -except dpath ; close all ;

% Set real theta1, optimal theta2, discrepancy version, whether modular
true_theta1 = 2;
discrep = 6;
if discrep < 5; opt_theta2=4/3 ; else ; opt_theta2=1 ; end
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs
xmin = 0.5 ; xrange = 0.5;
t2min = 0;
t2range = 5;

% Load previously gathered results
discrep_str = int2str(discrep);
if discrep == 7 ;  discrep_str = ['5_inf'] ; discrep = 5 ; end
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',discrep_str]);
load(locstr);

% Gather post-burn-in results
burn_in = results.settings.burn_in+1;
theta1 = results.theta1(burn_in:end,:);
theta2 = results.theta2(burn_in:end,:);
obs_rho = results.obs_rho(burn_in:end,:);
obs_lambda = results.obs_lambda(burn_in:end,:);
des_rho = results.des_rho(burn_in:end,:);
des_lambda = results.des_lambda(burn_in:end,:);

% Recover observations and target outcomes
obs_x  = results.settings.obs_x;
obs_y  = results.settings.obs_y;
obs_t2 = results.settings.obs_t2;
des_y  = results.settings.des_y;
m = size(theta1,1);
mean_y = results.settings.mean_y;
std_y = results.settings.std_y;

% Set desired observations (chosen for Gaussian quadrature)
n = 10 ; % Number of points to use for integration
[x,w] = lgwt(n,0,1); % Get points and weights


% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-5*eye(size(X)); % adds nugget for computational stablty
prior_mean = results.settings.mean_obs;
prior_cov = @(rho,x,t2,xp,t2p,lambda) ... 
    gp_cov(rho,[x ones(size(x,1),1).*t2],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,x,t2,xnew,t2new,rho,lambda) ...
    prior_mean(xnew,t2new) + ...
    (add_nug(prior_cov(rho,xnew,t2new,x,t2,lambda)) / ...
    add_nug(prior_cov(rho,x,t2,x,t2,lambda))) * ...
    (y - prior_mean(x,t2)) ;
updated_cov = @(x,t2,xnew,t2new,rho,lambda) ...
    add_nug(prior_cov(rho,xnew,t2new,xnew,t2new,lambda)) - ...
    (add_nug(prior_cov(rho,xnew,t2new,x,t2,lambda)) / ...
    add_nug(prior_cov(rho,x,t2,x,t2,lambda))) * ...
    add_nug(prior_cov(rho,x,t2,xnew,t2new,lambda)) ;
% clear results ; % No longer needed

% Get computer model output for each draw from the posterior (at des_x),
% and also get true output
comp_model_output = dual_calib_example_fn(repmat(x,m,1),xmin,xrange,...
    repelem(theta1,n,1),0,1,repelem(theta2,n,1),0,1,0,1,0);
true_output = dual_calib_example_fn(x,xmin,xrange,...
    true_theta1,0,1,opt_theta2,0,1,0,1,discrep);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
comp_model_output = reshape(comp_model_output,n,m)';
true_output = true_output';

% Add posterior mean of observation discrepancy
discrep_gp_post_means_std = comp_model_output * 0 ; % Pre-allocate space
discrep_gp_post_sds_std = comp_model_output * 0 ; % Pre-allocate space
obs_disc_std = zeros(m,size(obs_x,1)) ; % Pre-allocate space
for idx = 1:m
    rho = obs_rho(idx,:) ; lambda = obs_lambda(idx,:) ; 
    t1 = theta1(idx,:); t2 = (theta2(idx,:)-t2min)./t2range;
    % obs_disc_std holds the discrepancy between the observed values of y
    % and the computer model output for each draw of theta_1, on the
    % standardized scale.
    obs_disc_std(idx,:)=obs_y - dual_calib_example_fn(obs_x,xmin,xrange,...
        t1,0,1,obs_t2,t2min,t2range,...
        mean_y,std_y,0);
    % Gather discrepancy mean on standardized scale:
    d =updated_mean(obs_disc_std(idx,:)',obs_x,obs_t2,x,t2,rho,lambda);
    discrep_gp_post_means_std(idx,:) = d ;
    % Now get standard deviations of d:
    d_std_cov = updated_cov(obs_x,obs_t2,x,t2,rho,lambda) ; 
    d_std_sds = sqrt(diag(d_std_cov)+0.05) ; 
    discrep_gp_post_sds_std(idx,:) = d_std_sds ; 
    if mod(idx,100)==0 ; disp(m-idx) ; end
end

% Take a look at the observed discrepancies and the GP discrepancy means
figure();
obs_disc_mean = mean(obs_disc_std);
scatter3(obs_x,obs_t2,obs_disc_mean)
hold on;
scatter3(...
    x,...
    ones(size(x))*mean((theta2-t2min)./t2range),...
    mean(discrep_gp_post_means_std))

% Add discrepancy means to computer model output, rescale sds
discrep_gp_post_sds = discrep_gp_post_sds_std * std_y ;% Put on orig scale
discrep_gp_post_means = discrep_gp_post_means_std * std_y;% Put on orig scl
posterior_preds = comp_model_output + discrep_gp_post_means;
posterior_lb = posterior_preds - 2*discrep_gp_post_sds; % lower bound
posterior_ub = posterior_preds + 2*discrep_gp_post_sds; % upper bound

% Get average across x, using Gaussian quadrature
figure();
posterior_preds_avg = posterior_preds * w;
posterior_lb_avg = posterior_lb * w;
posterior_ub_avg = posterior_ub * w;
true_optimum = true_output * w;

% Make histograms
histogram(posterior_preds_avg); hold on;
histogram(posterior_lb_avg); histogram(posterior_ub_avg); 
ylims = get(gca,'Ylim');
plot([1 1]*true_optimum,ylims,'LineWidth',2);

% Show posterior predictive distribution of average output with true optim
posterior_distro = ...
    normpdf(linspace(0,1),posterior_preds,discrep_gp_post_sds);
