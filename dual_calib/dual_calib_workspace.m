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
discrep = 4;
modular = false;
obs_discrep = false; % Whether or not to include discrep term for real obs

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
avg_disc=0 ; % Actually nevermind
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
    'mean_y',mean_y,'std_y',std_y,'M',4e3,'burn_in',.5,...
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
    '2019-06-21_DCTO_obs-disc-prior0_discrep',int2str(discrep)]);
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
    '2019-06-21_KOH_calib_obs-disc-prior0_discrep',int2str(discrep)]);
% save(locstr,'results');

%% Perform CTO after KOH calibration, compare to DCTO
clc ; clearvars -except dpath ; close all ; 

% Which discrepancy version are we working with:
discrep = 4;

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
    '2019-06-21_KOH_calib_obs-disc-prior0_discrep',int2str(discrep)]);
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
    'mean_y',mean_y,'std_y',std_y,'M',1e4,'burn_in',.4,...
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
    '2019-06-21_DCTO_obs-disc-prior0_discrep',int2str(discrep)]);
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
results = KOH_results;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-06-21_KOH_calib_step2_obs-disc-prior0_discrep',...
    int2str(discrep)]);
% save(locstr,'results');

%% Find vars and MSEs and ADs of DCTO and KOH for each of 7 versions
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
t1_ad  = zeros(2,8);
t2_ad  = zeros(2,8);


% Loop through each discrep value and find MSEs
for ii=0:6%7
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
        '2019-06-21_DCTO_obs-disc-prior0_discrep',discrep_str]);
    koh1str = sprintf(['C:\\Users\\carle\\Documents',...
        '\\MATLAB\\NSF DEMS\\Phase 1\\',...
        'dual_calib\\dual_calib_stored_data\\'...
        '2019-06-21_KOH_calib_obs-disc-prior0_discrep',discrep_str]);
    koh2str = sprintf(['C:\\Users\\carle\\Documents',...
        '\\MATLAB\\NSF DEMS\\Phase 1\\',...
        'dual_calib\\dual_calib_stored_data\\'...
        '2019-06-21_KOH_calib_step2_obs-disc-prior0_discrep',discrep_str]);
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
    koh_t1_mse  = sum((koh_t1 - theta1).^2)/n ;
    koh_t2_mse  = sum((koh_t2 - theta2).^2)/n ;
    t1_mse(:,ii+1) = [dcto_t1_mse koh_t1_mse] ;
    t2_mse(:,ii+1) = [dcto_t2_mse koh_t2_mse] ;
    
    % Also var
    dcto_t1_var = var(dcto_t1) ; 
    dcto_t2_var = var(dcto_t2) ;    
    koh_t1_var  = var(koh_t1)  ;
    koh_t2_var  = var(koh_t2)  ;
    t1_var(:,ii+1) = [dcto_t1_var koh_t1_var] ;
    t2_var(:,ii+1) = [dcto_t2_var koh_t2_var] ;
    
    % And absolute deviation
    dcto_t1_ad = abs(mean(dcto_t1) - theta1);
    dcto_t2_ad = abs(mean(dcto_t2) - theta2);
    koh_t1_ad  = abs(mean(koh_t1) - theta1) ;
    koh_t2_ad  = abs(mean(koh_t2) - theta2) ;
    t1_ad(:,ii+1) = [dcto_t1_ad koh_t1_ad]  ;
    t2_ad(:,ii+1) = [dcto_t2_ad koh_t2_ad]  ;
    
end

t1_mse
t2_mse
t1_var
t2_var
t1_ad
t2_ad

%% Take a look at realizations of observation discrepancies
clc ; clearvars -except dpath ; close all ;

%%% Open figure and set which rho,lambda draws are used
f=figure('pos',[680 140 600 500]);
idx = randsample(4000,1);

%%% Set which discrepancy version is being examined
discrep = 4; 

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
% Unfinished?
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

% Set real theta1, optimal theta2, discrepancy version, whether modular,
% and known observation variance
sig2=0.05;
true_theta1 = 2;
discrep = 1;
if discrep < 5; opt_theta2=4/3 ; else ; opt_theta2=1 ; end
t1min = 1.5; t1range = 3 ;
t2min = 0  ; t2range = 5 ;
xmin  = .5 ; xrange  = .5;

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

% Recover observations and target outcomes
obs_x  = results.settings.obs_x;
obs_y  = results.settings.obs_y;
obs_t2 = results.settings.obs_t2;
des_y  = results.settings.des_y;
m = size(theta1,1);
mean_y = results.settings.mean_y;
std_y = results.settings.std_y;

% Set desired observations
n = 8 ; % Number of points to use for integration
[x,w] = lgwt(n,0,1); % Get points and weights


% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-5*eye(size(X)); % adds nugget for computational stablty
sig2=0.05; % known observation variance
prior_mean = results.settings.mean_obs;
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    prior_mean(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
    add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda)+sig2*eye(size(xo,1))))*...
    (y - prior_mean(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
    add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda)+sig2*eye(size(xo,1))))*...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
% clear results ; % No longer needed

% Get computer model output for each draw from the posterior (at x),
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
    d = updated_mean(obs_disc_std(idx,:)',obs_x,obs_t2,x,t2,rho,lambda);
    discrep_gp_post_means_std(idx,:) = d ;
    % Now get standard deviations of d:
    d_std_cov = updated_cov(obs_x,obs_t2,x,t2,rho,lambda) ; 
    d_std_sds = sqrt(diag(d_std_cov)+0.0) ; % Could add sig2 here
    discrep_gp_post_sds_std(idx,:) = d_std_sds ; 
    if mod(idx,1000)==0 ; disp(m-idx) ; end
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

% Generate samples of theta1,theta2 from prior distribution
theta1 = rand(m,1) * t1range + t1min;
theta2 = rand(m,1) * t2range + t2min;

% Get predicted output at each x point
prior_model_output = dual_calib_example_fn(repmat(x,m,1),xmin,xrange,...
    repelem(theta1,n,1),0,1,repelem(theta2,n,1),0,1,0,1,0);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
prior_model_output = reshape(prior_model_output,n,m)';

% Show posterior predictive distribution of average output with true optim,
% and include prior predictive distribution
figure();
for ii = 1:n
    posterior_distro = ...
        mean(normpdf(linspace(0,1),...
        posterior_preds(:,ii),discrep_gp_post_sds(:,ii)));
    subplot(2,n/2,ii);
    histogram(prior_model_output(:,ii),'Normalization','pdf',...
        'DisplayStyle','stairs'); hold on;
    plot(linspace(0,1),posterior_distro); 
    ylims = get(gca,'YLim');
    plot(true_output(ii)*[1 1],ylims,'LineWidth',2);
end


%% Get prior predictive distribution
clc ; clearvars -except dpath ; close all ;

% Set theta1, theta2, x ranges and mins
t1min = 1.5; t1range = 3 ;
t2min = 0  ; t2range = 5 ;
xmin  = .5 ; xrange  = .5;

% Number of samples and number of x points at which to get output:
m = 16000; % samples
n = 10 ; % x points
x = linspace(0,1,n)';

% Generate samples of theta1,theta2 from prior distribution
theta1 = rand(m,1) * t1range + t1min;
theta2 = rand(m,1) * t2range + t2min;

% Get predicted output at each x point
prior_model_output = dual_calib_example_fn(repmat(x,m,1),xmin,xrange,...
    repelem(theta1,n,1),0,1,repelem(theta2,n,1),0,1,0,1,0);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
prior_model_output = reshape(prior_model_output,n,m)';

% Show posterior predictive distribution of average output with true optim
figure();
for ii = 1:n
    subplot(2,n/2,ii);
    histogram(prior_model_output(:,ii)); 
end

%% Get predictive distribution from KOH+CTO/DCTO, compare to Pareto front
clc ; clearvars -except dpath ; close all ;

% Set real theta1, optimal theta2, discrepancy version, whether modular,
% and known observation variance
sig2=0.05;
true_theta1 = 2;
discrep = 2;
if discrep < 5; opt_theta2=4/3 ; else ; opt_theta2=1 ; end
t1min = 1.5; t1range = 3 ;
t2min = 0  ; t2range = 5 ;
xmin  = .5 ; xrange  = .5;

% Load previously gathered results
discrep_str = int2str(discrep);
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
burn_in=results.settings.burn_in; 
% Get theta1, theta2, rho, lambda, observations and targets from dcto
dcto_t1 = results.theta1(burn_in+1:end); n = length(dcto_t1);
dcto_t2 = results.theta2(burn_in+1:end);
dcto_rho = results.obs_rho(burn_in+1:end,:);
dcto_lambda = results.obs_lambda(burn_in+1:end,:);
% Recover observations and target outcomes
obs_x  = results.settings.obs_x;
obs_y  = results.settings.obs_y;
obs_t2 = results.settings.obs_t2;
des_y  = results.settings.des_y;
m = size(dcto_t1,1);
mean_y = results.settings.mean_y;
std_y = results.settings.std_y;
prior_mean = results.settings.mean_obs;
% Get theta1, theta2, rho, lambda from koh and cto
load(koh1str); koh_t1 = randsample(results.theta1(burn_in+1:end),m);
koh_rho = results.obs_rho(burn_in+1:end,:);
koh_lambda = results.obs_lambda(burn_in+1:end,:);
load(koh2str); koh_t2 = randsample(results.theta1(burn_in+1:end),m);


% Set desired observations
n = 1 ; % Number of points to use for integration
x = linspace(0,1,n)'; % Get points and weights
% [x,w] = lgwt(n,0,1); % Get points and weights
 

% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-5*eye(size(X)); % adds nugget for computational stablty
sig2=0.05; % known observation variance
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    prior_mean(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
    add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda)+sig2*eye(size(xo,1))))*...
    (y - prior_mean(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
    add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda)+sig2*eye(size(xo,1))))*...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
% clear results ; % No longer needed

% Get computer model output for each draw from the posterior (at x),
% and also get true output
comp_model_output_dcto =dual_calib_example_fn(repmat(x,m,1),xmin,xrange,...
    repelem(dcto_t1,n,1),0,1,repelem(dcto_t2,n,1),0,1,0,1,0);
comp_model_output_koh =dual_calib_example_fn(repmat(x,m,1),xmin,xrange,...
    repelem(koh_t1,n,1),0,1,repelem(koh_t2,n,1),0,1,0,1,0);
true_output = dual_calib_example_fn(x,xmin,xrange,...
    true_theta1,0,1,opt_theta2,0,1,0,1,discrep);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
comp_model_output_dcto = reshape(comp_model_output_dcto,n,m)';
comp_model_output_koh = reshape(comp_model_output_koh,n,m)';
true_output = true_output';

% Add posterior mean of observation discrepancy
discrep_gp_post_means_std_dcto = ...
    comp_model_output_dcto * 0 ;%Pre-allocate space
discrep_gp_post_sds_std_dcto = ...
    comp_model_output_dcto * 0 ; % Pre-allocate space
obs_disc_std_dcto = zeros(m,size(obs_x,1)) ; % Pre-allocate space
discrep_gp_post_means_std_koh = ...
    comp_model_output_koh * 0 ; %Pre-allocate space
discrep_gp_post_sds_std_koh = ...
    comp_model_output_koh * 0 ; % Pre-allocate space
obs_disc_std_koh = zeros(m,size(obs_x,1)) ; % Pre-allocate space
for idx = 1:m
    rho_dcto = dcto_rho(idx,:) ; lambda_dcto = dcto_lambda(idx,:) ; 
    t1_dcto = dcto_t1(idx,:); t2_dcto = (dcto_t2(idx,:)-t2min)./t2range;
    % obs_disc_std holds the discrepancy between the observed values of y
    % and the computer model output for each draw of theta_1, on the
    % standardized scale.
    obs_disc_std_dcto(idx,:)= ...
        obs_y - dual_calib_example_fn(obs_x,xmin,xrange,...
        t1_dcto,0,1,obs_t2,t2min,t2range,...
        mean_y,std_y,0);
    % Gather discrepancy mean on standardized scale:
    d_dcto=updated_mean(obs_disc_std_dcto(idx,:)',obs_x,obs_t2,...
        x,t2_dcto,rho_dcto,lambda_dcto);
    discrep_gp_post_means_std_dcto(idx,:) = d_dcto ;
    % Now get standard deviations of d:
    d_std_cov_dcto = ...
        updated_cov(obs_x,obs_t2,x,t2_dcto,rho_dcto,lambda_dcto) ; 
    d_std_sds_dcto = sqrt(diag(d_std_cov_dcto)+0.0) ; % Could add sig2 here
    discrep_gp_post_sds_std_dcto(idx,:) = d_std_sds_dcto ; 
    % Now do it all for KOH
    rho_koh = koh_rho(idx,:) ; lambda_koh = koh_lambda(idx,:) ; 
    t1_koh = koh_t1(idx,:); t2_koh = (koh_t2(idx,:)-t2min)./t2range;
    % obs_disc_std holds the discrepancy between the observed values of y
    % and the computer model output for each draw of theta_1, on the
    % standardized scale.
    obs_disc_std_koh(idx,:)= ...
        obs_y - dual_calib_example_fn(obs_x,xmin,xrange,...
        t1_koh,0,1,obs_t2,t2min,t2range,...
        mean_y,std_y,0);
    % Gather discrepancy mean on standardized scale:
    d_koh=updated_mean(obs_disc_std_koh(idx,:)',obs_x,obs_t2,...
        x,t2_koh,rho_koh,lambda_koh);
    discrep_gp_post_means_std_koh(idx,:) = d_koh ;
    % Now get standard deviations of d:
    d_std_cov_koh = ...
        updated_cov(obs_x,obs_t2,x,t2_koh,rho_koh,lambda_koh) ; 
    d_std_sds_koh = sqrt(diag(d_std_cov_koh)+0.0) ; % Could add sig2 here
    discrep_gp_post_sds_std_koh(idx,:) = d_std_sds_koh ; 
    if mod(idx,1000)==0 ; disp(m-idx) ; end
end

% Take a look at the observed discrepancies and the GP discrepancy means
figure();
obs_disc_mean_dcto = mean(obs_disc_std_dcto);
scatter3(obs_x,obs_t2,obs_disc_mean_dcto)
hold on;
scatter3(...
    x,...
    ones(size(x))*mean((dcto_t2-t2min)./t2range),...
    mean(discrep_gp_post_means_std_dcto))
obs_disc_mean_koh = mean(obs_disc_std_koh);
scatter3(obs_x,obs_t2,obs_disc_mean_koh)
scatter3(...
    x,...
    ones(size(x))*mean((koh_t2-t2min)./t2range),...
    mean(discrep_gp_post_means_std_koh))

% Add discrepancy means to computer model output, rescale sds
discrep_gp_post_sds_dcto = ...
    discrep_gp_post_sds_std_dcto * std_y ;% Put on orig scale
discrep_gp_post_means_dcto = ...
    discrep_gp_post_means_std_dcto * std_y;% Put on orig scl
posterior_preds_dcto = comp_model_output_dcto + discrep_gp_post_means_dcto;
posterior_lb_dcto = ...
    posterior_preds_dcto - 2*discrep_gp_post_sds_dcto; % lower bound
posterior_ub_dcto = ...
    posterior_preds_dcto + 2*discrep_gp_post_sds_dcto; % upper bound
discrep_gp_post_sds_koh = ...
    discrep_gp_post_sds_std_koh * std_y ;% Put on orig scale
discrep_gp_post_means_koh = ...
    discrep_gp_post_means_std_koh * std_y;% Put on orig scl
posterior_preds_koh = comp_model_output_koh + discrep_gp_post_means_koh;
posterior_lb_koh = ...
    posterior_preds_koh - 2*discrep_gp_post_sds_koh; % lower bound
posterior_ub_koh = ...
    posterior_preds_koh + 2*discrep_gp_post_sds_koh; % upper bound

% Get average across x, using Gaussian quadrature
% figure();
% posterior_preds_avg_dcto = posterior_preds_dcto * w;
% posterior_lb_avg_dcto = posterior_lb_dcto * w;
% posterior_ub_avg_dcto = posterior_ub_dcto * w;
% posterior_preds_avg_koh = posterior_preds_koh * w;
% posterior_lb_avg_koh = posterior_lb_koh * w;
% posterior_ub_avg_koh = posterior_ub_koh * w;
% true_optimum = true_output * w;
% % Make histograms
% histogram(posterior_preds_avg_dcto); hold on;
% histogram(posterior_lb_avg_dcto); histogram(posterior_ub_avg_dcto); 
% histogram(posterior_preds_avg_koh);
% histogram(posterior_lb_avg_koh); histogram(posterior_ub_avg_koh); 
% ylims = get(gca,'Ylim');
% plot([1 1]*true_optimum,ylims,'LineWidth',2);

% Generate samples of theta1,theta2 from prior distribution
theta1 = rand(m,1) * t1range + t1min;
theta2 = rand(m,1) * t2range + t2min;

% Get predicted output at each x point
prior_model_output = dual_calib_example_fn(repmat(x,m,1),xmin,xrange,...
    repelem(theta1,n,1),0,1,repelem(theta2,n,1),0,1,0,1,0);

% Reshape so that each row corresponds to a single draw of (theta1,theta2)
prior_model_output = reshape(prior_model_output,n,m)';

% Show posterior predictive distribution of average output with true optim,
% and include prior predictive distribution
figure();
xpred = linspace(0,2.5,200); % Support of posterior pred distr
for ii = 1:n
    posterior_distro_dcto = ...
        mean(normpdf(xpred,...
        posterior_preds_dcto(:,ii),discrep_gp_post_sds_dcto(:,ii)));
    posterior_distro_koh = ...
        mean(normpdf(xpred,...
        posterior_preds_koh(:,ii),discrep_gp_post_sds_koh(:,ii)));
    subplot(1,n,ii);
    histogram(prior_model_output(:,ii),'Normalization','pdf',...
        'DisplayStyle','stairs'); hold on;
    plot(xpred,posterior_distro_dcto);
    plot(xpred,posterior_distro_koh); 
    ylims = get(gca,'YLim');
    plot(true_output(ii)*[1 1],ylims,'LineWidth',2);
    histogram(posterior_preds_dcto(:,ii),'Normalization','pdf',...
        'DisplayStyle','stairs')
    
    % Get posterior variances
    posterior_mean_dcto = xpred * posterior_distro_dcto' / ...
        sum(posterior_distro_dcto);
    posterior_var_dcto = sum(xpred.^2 .* ... 
        posterior_distro_dcto)/sum(posterior_distro_dcto) - ...
        posterior_mean_dcto^2;
    posterior_mean_koh = xpred * posterior_distro_koh' / ...
        sum(posterior_distro_koh);
    posterior_var_koh = sum(xpred.^2 .* ... 
        posterior_distro_koh)/sum(posterior_distro_koh) - ...
        posterior_mean_koh^2;
end

%% Compare DCTO/KOH+CTO on g2 at specified setting for (c,d)
% To investigate the unusually big difference between DCTO and KOH+CTO for
% the case of g2 with (c,d) = (.65, .075)
clc ; clearvars -except dpath ; close all ; 

% Set c,d
c = 0.65;
d = 0.075;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, discrepancy version, whether modular
theta1 = 2;
discrep = 4;
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs

% Set des_x size
des_x_size = 50;

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
% Only use specified c,d if we're using discrep 4
if discrep == 4
    obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
        xmin,xrange,...
        (obs_t1-t1min)/t1range,t1min,t1range,...
        (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep,c,d);
else
    obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
        xmin,xrange,...
        (obs_t1-t1min)/t1range,t1min,t1range,...
        (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);
end

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
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
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 
% The below lines are included so that an informative discrepancy mean can
% be set for discrep=3. 
% c = 0.055 ; % This controls the size of the additive discrepancy
% d = 0; %Raising/lowering c makes discrepancy more/less aggressive
% obs_discrep_mean = @(xx,tt2) ...
%     (c*(xx*xrange+xmin).*(tt2*t2range+t2min)+d)/std_y;

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform dual calibration
DCTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform KOH calibration

% Get "real" observations
% Get same observations used in DCTO version:
obs_y = DCTO_results.settings.obs_y * std_y + mean_y;

% Our x is now bivariate, incorporating the former x plus now also what was
% t2. This is because t2 is treated here as a control variable, not
% something over which we are attempting to optimize.

% Now set desired observations; since we're doing KOH, there are none
des_x = zeros(0,2); des_y = [];

% Define objective function
mean_sim = @(a,b,c) dual_calib_example_fn(...
    a(:,1),xmin,xrange,...
    b,t1min,t1range,... 
    a(:,2),t2min,t2range,...
    mean_y,std_y);


% Get settings
settings = MCMC_dual_calib_settings(zeros(0,2),sim_t1,[],sim_y,...
    [obs_x,obs_t2],[],obs_y,des_x,des_y,'min_x',[xmin,t2min],...
    'range_x',[xrange,t2range],...
    'min_t1',t1min,'range_t1',t1range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',modular,...
    'EmulatorMean',mean_sim);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
KOH_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform CTO after KOH

% Get theta1 dist from KOH:
theta1_samps = ...
    (KOH_results.theta1(KOH_results.settings.burn_in:end) - ...
    KOH_results.settings.min_t1) / KOH_results.settings.range_t1 ;

% Define objective function
mean_sim = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,...
    randsample(theta1_samps,1),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)
%TEMP: give set t1 value 
% mean_sim = @(a,b,c) dual_calib_example_fn(...
%     a,xmin,xrange,1/12,t1min,t1range,...
%     b,t2min,t2range,...
%     mean_y,std_y);

% Now set desired observations
obs_x = linspace(0,1,des_x_size)' * xrange + xmin;
obs_y = zeros(size(obs_x,1),1);
obs_t2 = [] ; % We set obs_t2 to be empty, since we only have a 
% KOH calibration here, no DCTO secondary calibration. Similarly we will
% set des_x,des_y to be empty
des_x = [] ; des_y = [];

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',true,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

% Get all results on original scale


% First, get prior and posterior theta1
f1 = figure('pos',[10 10 550 225]);
lcol = [218 165 32]/255 ; % Color for line
subplot(1,2,1);
% Plot prior
fill([t1min t1min + t1range t1min + t1range t1min],...
    [0 0 1/t1range 1/t1range],'g','EdgeColor','none');
xlim([t1min t1min + t1range]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
% Save it:
set(f1,'color','white');
% saveas(f1,'FIG_dual_calib_post_theta1-2.png');

% Second, get prior and posterior theta2
% f2 = figure('pos',[10 320 400 300]);
subplot(1,2,2);
% left_color = [.5 .5 0];
% right_color = [.5 0 .5];
% set(f2,'defaultAxesColorOrder',[left_color; right_color]);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
if discrep == 4
    fmfn =@(z) dual_calib_example_fn(...
        .75,0,1,theta1,0,1,z,0,1,0,1,discrep,c,d);
else
    fmfn =@(z) dual_calib_example_fn(...
        .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
end
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
% title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
% set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(5),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
% all_results_locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-07-19_KOHCTO_DCTO_comparisons']);
% load(all_results_locstr);
% idx = max(size(all_results)) + 1;
% all_results{idx} = struct('cd',[c,d],'discrep',discrep,...
%     'KOH_results',KOH_results,'CTO_results',CTO_results,...
%     'DCTO_results',DCTO_results);
% save(all_results_locstr,'all_results');
discrep_str = int2str(discrep);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-01_DCTO_discrep',discrep_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-01_KOH_discrep_',discrep_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-01_CTO_after_KOH_discrep_',discrep_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Compare DCTO and KOH+CTO with calib of t2 and design on t1
clc ; clearvars -except dpath ; close all ; 

% Set c,d
c = 0.65;
d = 0.075;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta2, discrepancy version, whether modular
theta2 = 4/3;
discrep = 4;
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs

% Set des_x size
des_x_size = 8;

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
obs_x = design.obs_x; 
% The saved design treats theta1 as the calibration input, and so contains
% observations at known theta2 levels. Since we're now treating theta2 as
% the calibration parameter, we'll convert this design so that it treats
% theta1 as known.
obs_t1 = ((design.obs_t2 - t2min)./t2range)*t1range + t1min; 
mean_y = design.mean_y;
std_y = design.std_y;

% Make a col vector based on true theta2
obs_t2 = ones(size(obs_x,1),1) * theta2;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep,c,d);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 
% The below lines are included so that an informative discrepancy mean can
% be set for discrep=3. 
% c = 0.055 ; % This controls the size of the additive discrepancy
% d = 0; %Raising/lowering c makes discrepancy more/less aggressive
% obs_discrep_mean = @(xx,tt2) ...
%     (c*(xx*xrange+xmin).*(tt2*t2range+t2min)+d)/std_y;

% Set objective function
mean_sim = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,c,t1min,t1range,...
    b,t2min,t2range,...
    mean_y,std_y);

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t1,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',t1min,'range_t2',t1range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular,...
    'EmulatorMean',mean_sim);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform dual calibration
DCTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform KOH calibration

% Get "real" observations
% Get same observations used in DCTO version:
obs_y = DCTO_results.settings.obs_y * std_y + mean_y;

% Now set desired observations; since we're doing KOH, there are none
des_x = []; des_y = [];

mean_sim = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,...
    c,t1min,t1range,...
    b,t2min,t2range,...
    mean_y,std_y);

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t1,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',t1min,'range_t2',t1range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',modular,...
    'EmulatorMean',mean_sim);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
KOH_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform CTO after KOH

% Get theta1 dist from KOH:
theta1_samps = ...
    (KOH_results.theta1(KOH_results.settings.burn_in:end) - ...
    KOH_results.settings.min_t1) / KOH_results.settings.range_t1 ;

% Define objective function
mean_sim = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,b,t1min,t1range,...
    randsample(theta1_samps,1),t2min,t2range,...
    mean_y,std_y);
%TEMP: give set t1 value 
% mean_sim = @(a,b,c) dual_calib_example_fn(...
%     a,xmin,xrange,1/12,t1min,t1range,...
%     b,t2min,t2range,...
%     mean_y,std_y);

% Now set desired observations
obs_x = linspace(0,1,des_x_size)' * xrange + xmin;
obs_y = zeros(size(obs_x,1),1);
obs_t2 = linspace(t1min, t1min+t1range,des_x_size)';
% obs_t2 = obs_x ; % obs_t2 doesn't matter, it is ignored by the obj fn
% Notice we named the above obs_x,obs_y rather than des_x,des_y. This is
% because we are performing KOH calibration on theta2. We set des_x,des_y
% to be empty, since we only have a KOH calibration here, no DCTO secondary
% calibration.

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t1min,'range_t2',t1range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',true,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim);
%'obs_discrep_mean',@(x,t)0.055*x.*t,...% Store this here

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

% Get all results on original scale


% First, get prior and posterior theta1
%f1 = figure('pos',[10 10 400 300]); 
f1 = figure('pos',[10 10 600 250]);
subplot(1,2,1);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.05);
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.05);
% Plot true theta2
plot([theta2 theta2],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it:
set(f1,'color','white');
% saveas(f1,'FIG_dual_calib_post_theta1-2.png');

% Second, get prior and posterior theta2
%f2 = figure('pos',[10 320 400 300]);
subplot(1,2,2)
left_color = [.5 .5 0];
right_color = [.5 0 .5];
% set(f2,'defaultAxesColorOrder',[left_color; right_color]);
% Plot prior
fill([t1min t1min + t1range t1min + t1range t1min],...
    [0 0 1/t1range 1/t1range],'g','EdgeColor','none');
xlim([t1min t1min + t1range]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(.75,0,1,z,0,1,theta2,0,1,0,1,discrep,c,d);
theta1 = fmincon(fmfn,2,[],[],[],[],t1min,t1min+t1range);
yyaxis left ;
plot([theta1 theta1],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','KOH','DCTO','Optimal value');
title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
% Save it as is:
set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');


%% Take compare DCTO target discrepancy with CTO discrepancy
clc ; clearvars -except dpath ; close all ;

%%% Open figure and set which rho,lambda draws are used
f=figure('pos',[680 10 500 400]);
idx = randsample(4000,1);

%%% Set which discrepancy version is being examined
discrep = 3; 

%%% Define input rescaling settings
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

%%% Set true/optimal values of theta1, theta2
true_theta1 = 2;
true_theta2 = 4/3;

%%% Load the DCTO MCMC results
discrep_str = int2str(discrep);
if discrep == 7 discrep_str = ['5_inf'] ; discrep = 5 ; end
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-01_DCTO_discrep',discrep_str]);
load(locstr);
burnin = DCTO_results.settings.burn_in;
rhos = DCTO_results.des_rho(burnin:end,:) ;
lambdas = DCTO_results.des_lambda(burnin:end,:);
theta1s = DCTO_results.theta1(burnin:end,:);
theta2s = DCTO_results.theta2(burnin:end,:);

%%% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-4*eye(size(X)); % adds nugget for computational stablty
prior_mean = @(x) zeros(size(x)); 
prior_cov = @(rho,x,xp,lambda) ... 
    gp_cov(rho,x,xp,lambda,false);
updated_mean = @(y,x,xnew,rho,lambda) ...
    prior_mean(xnew) + ...
    add_nug(prior_cov(rho,xnew,x,lambda)) * ...
    inv(add_nug(prior_cov(rho,x,x,lambda))) * ...
    (y - prior_mean(x)) ;
updated_cov = @(x,xnew,rho,lambda) ...
    add_nug(prior_cov(rho,xnew,xnew,lambda)) - ...
    add_nug(prior_cov(rho,xnew,x,lambda)) * ...
    inv(add_nug(prior_cov(rho,x,x,lambda))) * ...
    add_nug(prior_cov(rho,x,xnew,lambda)) ;

%%% Get vectors for observed points and prediction points
des_x = DCTO_results.settings.des_x ; 
des_y = DCTO_results.settings.des_y ; 
h=41; % Grid mesh setting for x
xp = linspace(0,1,h);
xp=xp(:) ;

%%% Get mean of discrep for a random draw of rho,lambda
rho = rhos(idx,:) ; lambda = lambdas(idx,:) ; 
theta1 = mean(theta1s) ; 
theta2 = mean(theta2s) ; % theta2s(idx,:) ;
[theta1 theta2] % Have a look
% desc_disc is the discrepancy we observe -- the difference between the
% desired observations and the computer model run using estimates of
% theta1, theta2. It is on the standardized scale.
des_disc = des_y - dual_calib_example_fn(des_x,xmin,xrange,...
    theta1,0,1,theta2,0,1,...
    DCTO_results.settings.mean_y,...
    DCTO_results.settings.std_y,0);
% Gather discrepancy mean on standardized scale:
d_std = updated_mean(des_disc,des_x,xp,rho,lambda) ;
% Transform to original scale:
d = d_std * DCTO_results.settings.std_y;
% Now get standard deviations of d:
d_std_cov = updated_cov(des_x,xp,rho,lambda) ; 
d_std_sds = sqrt(diag(d_std_cov)+0.05) ; 
d_sds = d_std_sds * DCTO_results.settings.std_y;

%%% Plot the discrepancy
Y = dual_calib_example_fn(xp,xmin,xrange,true_theta1,0,1,...
    true_theta2,0,1,0,1,0);
Yd= zeros(length(xp),1);

% Take a look
% Discrep is the true discrepancy between 1. the computer model with the 
% true value of theta1 and the optimal value of theta2, and 2. the desired
% observations (in this case, constant 0).
Discrep = Yd-Y; 
plot(xp,Discrep,'b','LineWidth',2);
hold on;
% Now we add the estimated discrepancy:
d_upper = d+2*d_sds ; d_lower = d-2*d_sds ;
plot(xp,d,'r','Linewidth',2);
plot(xp,d_upper,'r--','LineWidth',2);
plot(xp,d_lower,'r--','LineWidth',2);

%%%%%%%%%%%%%%%%%%%
%%% Now do the same for KOH+CTO
idx = randsample(4000,1);

%%% Load the KOH MCMC results
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-01_KOH_discrep_',discrep_str]);
load(locstr);
theta1s = KOH_results.theta1(burnin:end,:);

%%% Load the CTO MCMC results
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-01_CTO_after_KOH_discrep_',discrep_str]);
load(locstr);
rhos = CTO_results.obs_rho(burnin:end,:) ;
lambdas = CTO_results.obs_lambda(burnin:end,:);
theta2s = CTO_results.theta1(burnin:end,:);

%%% Get vectors for observed points and prediction points
des_x = CTO_results.settings.obs_x ; 
des_y = CTO_results.settings.obs_y ; 

%%% Get mean of discrep for a random draw of rho,lambda
rho = rhos(idx,:) ; lambda = lambdas(idx,:) ; 
theta1 = mean(theta1s) ; 
theta2 = mean(theta2s) ; %theta2s(idx,:) ;
% desc_disc is the discrepancy we observe -- the difference between the
% desired observations and the computer model run using estimates of
% theta1, theta2. It is on the standardized scale.
des_disc = des_y - dual_calib_example_fn(des_x,xmin,xrange,...
    theta1,0,1,theta2,0,1,...
    CTO_results.settings.mean_y,CTO_results.settings.std_y,0);
% Gather discrepancy mean on standardized scale:
d_std = updated_mean(des_disc,des_x,xp,rho,lambda) ;
% Transform to original scale:
d = d_std * CTO_results.settings.std_y ;
% Now get standard deviations of d:
d_std_cov = updated_cov(des_x,xp,rho,lambda) ; 
d_std_sds = sqrt(diag(d_std_cov)+0.05) ; 
d_sds = d_std_sds * CTO_results.settings.std_y;

% Now we add the estimated discrepancy to the plot:
d_upper = d+2*d_sds ; d_lower = d-2*d_sds ;
plot(xp,d,'g','Linewidth',2);
plot(xp,d_upper,'g--','LineWidth',2);
plot(xp,d_lower,'g--','LineWidth',2);

[theta1 theta2]