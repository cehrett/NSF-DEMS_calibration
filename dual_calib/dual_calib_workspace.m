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
    T2(:),t2min,t2range,0,1,0,true),length(x),length(t1),length(t2));

% Take a look
xidx=100;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
surf(tt1 * t1range + t1min,tt2*t2range+t2min,reshape(Y(:,xidx,:),100,100));

%% Gather data to use for emulator and for "real" observations
clc ; clearvars -except dpath ; close all ;

% Set discrepancy version
discrep = 0; % If collecting raw data for emulator this should be 0

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Get simulation design points and observations
X = lhsdesign(100,3);
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;
% Get simulation observations
sim_y = dual_calib_example_fn(sim_x,sim_t1,sim_t2,discrep);

% % Now set "real" theta1 and get "real" observations
% theta1 = 2;
% X = lhsdesign(30,2);
% obs_x = X(:,1) * xrange + xmin;
% obs_t2 = X(:,2) * t2range + t2min;
% % Make a col vector based on true theta1
% obs_t1 = ones(size(obs_x,1),1) * theta1;
% % Get "real" observations without noise or discrepancy
% obs_y_noiseless = dual_calib_example_fn(obs_x,obs_t1,obs_t2);
% 
% % Now noise it up (still keep no discrepancy for now)
% sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% std_y = std(sim_y);
% mean_y = mean(sim_y);
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% 
% % Now set desired observations
% des_x = linspace(0,1,8)' * xrange + xmin;
% des_y = zeros(size(des_x,1),1);

% Now package everything up and save it
clear X sigma obs_y_noiseless dpath
save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-09-10_dual_calib_raw_data_for_emulator']);

% Get dpath back
direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

%% Get MLEs for covariance parameters for emulator GP
clc ; clearvars -except dpath ; close all ; 

% Set discrepancy
discrep = 1;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-08-28_dual_calib_raw_data_discrep',int2str(discrep)]);

% Normalize the simulation inputs, standardize the outputs
x = (sim_x - xmin) ./ xrange ;
t1 = (sim_t1 - t1min) ./ t1range;
t2 = (sim_t2 - t2min) ./ t2range;
y = (sim_y - mean(sim_y)) ./ std(sim_y);

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
inp


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
discrep = 1;
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


%% Compare DCTO target discrepancy with CTO discrepancy
clc ; clearvars -except dpath ; close all ;

%%% Open figure and set which rho,lambda draws are used
f=figure('pos',[680 10 500 400]);
idx = randsample(4000,1);

%%% Set which discrepancy version is being examined, which # of des_x
discrep = 2; 
des_x_size = 3;

%%% Set nugget size
nugsize = 0.000005;

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
des_x_size_str = int2str(des_x_size);
if discrep == 7 discrep_str = ['5_inf'] ; discrep = 5 ; end
% locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-08-01_DCTO_discrep',discrep_str]);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
load(DCTO_locstr);
burnin = DCTO_results.settings.burn_in;
rhos = DCTO_results.des_rho(burnin:end,:) ;
lambdas = DCTO_results.des_lambda(burnin:end,:);
theta1s = DCTO_results.theta1(burnin:end,:);
theta2s = DCTO_results.theta2(burnin:end,:);

%%% Define the updated mean and covariance functions for the discrepancy
add_nug = @(X) X+1e-10*eye(size(X)); % adds nugget for computatnl stablty
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
% h=length(des_x); % Grid mesh setting for x
h = 50;
xp = linspace(0,1,h);
xp=xp(:) ;

%%% Get mean of discrep for a random draw or mean of rho,lambda
% rho = rhos(idx,:) ; lambda = lambdas(idx,:) ; 
rho = mean(rhos) ; lambda = mean(lambdas);
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
d_std_sds = sqrt(diag(d_std_cov) + nugsize) ; 
d_sds = d_std_sds * DCTO_results.settings.std_y;

%%% Plot the true discrepancy
Y = dual_calib_example_fn(xp,xmin,xrange,true_theta1,0,1,...
    true_theta2,0,1,0,1,0);
Yd= zeros(length(xp),1);

% Take a look
% Now we add the estimated discrepancy:
d_upper = d+2*d_sds ; d_lower = d-2*d_sds ;
plot(xp,d,'r','Linewidth',2);
hold on;
plot(xp,d_upper,'r--','LineWidth',2);
plot(xp,d_lower,'r--','LineWidth',2);

%%%%%%%%%%%%%%%%%%%
%%% Now do the same for KOH+CTO
idx = randsample(4000,1);

%%% Load the KOH MCMC results
% locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-08-01_KOH_discrep_',discrep_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);

load(KOH_locstr);
theta1s = KOH_results.theta1(burnin:end,:);

%%% Load the CTO MCMC results
% locstr = sprintf(['C:\\Users\\carle\\Documents',...
%     '\\MATLAB\\NSF DEMS\\Phase 1\\',...
%     'dual_calib\\dual_calib_stored_data\\'...
%     '2019-08-01_CTO_after_KOH_discrep_',discrep_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
load(CTO_locstr);
rhos = CTO_results.obs_rho(burnin:end,:) ;
lambdas = CTO_results.obs_lambda(burnin:end,:);
theta2s = CTO_results.theta1(burnin:end,:);

%%% Get vectors for observed points and prediction points
des_x = CTO_results.settings.obs_x ; 
des_y = CTO_results.settings.obs_y ; 

%%% Get mean of discrep for a random draw or mean of rho,lambda
% rho = rhos(idx,:) ; lambda = lambdas(idx,:) ; 
rho = mean(rhos) ; lambda = mean(lambdas) ; 
theta1 = mean(theta1s) ; 
theta2 = mean(theta2s) ; %theta2s(idx,:) ;
% des_disc is the discrepancy we observe -- the difference between the
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
d_std_sds = sqrt(diag(d_std_cov) + nugsize) ; 
d_sds = d_std_sds * CTO_results.settings.std_y;

% Now we add the estimated discrepancy to the plot:
d_upper = d+2*d_sds ; d_lower = d-2*d_sds ;
plot(xp,d,'g','Linewidth',2);
plot(xp,d_upper,'g--','LineWidth',2);
plot(xp,d_lower,'g--','LineWidth',2);

% Now add the true discrepancy
% Discrep is the true discrepancy between 1. the computer model with the 
% true value of theta1 and the optimal value of theta2, and 2. the desired
% observations (in this case, constant 0).
Discrep = Yd-Y; 
plot(xp,Discrep,'b','LineWidth',2);

[theta1 theta2]

%% Compare posterior dists from DCTO, KOH+CTO when obs disc not used
clc ; clearvars -except dpath ; close all

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, discrepancy version, whether modular
theta1 = 2;
discrep = 0;
modular = false;
obs_discrep = false; % Whether or not to include discrep term for real obs

% Set des_x size
des_x_size = 20;

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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

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

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'emulator',false,...
    'modular',modular);

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
    'des_discrep',false,...
    'emulator',false,...
    'modular',modular,...
    'EmulatorMean',mean_sim);

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
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim);

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
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
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
suptitle('DCTO and KOH+CTO, modeled without observation discrepancy');
flushLegend(lg1,f1.Children(5),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_DCTO_discrep',discrep_str,'_no_obs_discr_GP']);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_KOH_discrep',discrep_str,'_no_obs_discr_GP']);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_CTO_after_KOH_discrep',discrep_str,'_no_obs_discr_GP']);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Compare DCTO and KOH+CTO with various true observation error variances
clc ; clearvars -except dpath ; close all ; 

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, discrepancy version, whether modular
theta1 = 2;
discrep = 4;
new_sigma = sqrt(0.2);
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs
known_var = true; % Whether or not to treat new_sigma as known
if known_var, sigma2=new_sigma^2; else sigma2=0.05 ; end

% Set des_x size
des_x_size = 20;

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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
% We have a saved set of noise, for consistency. This saved set of noise is
% based on the old_sigma value for variance of 0.05. We're now going to use
% this saved noise to test other levels of error variance. So we're going
% to load the same saved noise and convert it to have different variance
% levels than the original.
old_sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations

% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
tempnoise = tempnoise / old_sigma * new_sigma;
obs_y = obs_y_noiseless + tempnoise;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'emulator',false,...
    'modular',modular,...
    'ObsVar',sigma2);

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
    'des_discrep',false,...
    'emulator',false,...
    'modular',modular,...
    'EmulatorMean',mean_sim,...
    'ObsVar',sigma2);

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
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim);

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
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
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
suptitle(sprintf('DCTO and KOH+CTO using observation error var %g',...
    new_sigma^2));
flushLegend(lg1,f1.Children(5),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_DCTO_discrep',discrep_str,'_err_var%g'],...
    new_sigma^2*100);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_KOH_discrep',discrep_str,'_err_var%g'],...
    new_sigma^2*100);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_CTO_after_KOH_discrep',discrep_str,'_err_var%g'],...
    new_sigma^2*100);
if known_var 
    DCTO_locstr = strcat(DCTO_locstr,'_known') ; 
    KOH_locstr = strcat(KOH_locstr,'_known') ; 
    CTO_locstr = strcat(CTO_locstr,'_known') ; 
end
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Examine previous results of DCTO and KOH+CTO with various sigma2 vals
clc ; clearvars -except dpath ; close all ;

% Select which discrepancy and sigma2 we want to look at
discrep = 1;
new_sigma = sqrt(0.01);

% Define inputs mins and ranges, true values
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;
theta1=2;
known_var = true; % Tells whether obs var treated as known

% Load the results
discrep_str = int2str(discrep);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_DCTO_discrep',discrep_str,'_err_var%g'],...
    new_sigma^2*100);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_KOH_discrep',discrep_str,'_err_var%g'],...
    new_sigma^2*100);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-02_CTO_after_KOH_discrep',discrep_str,'_err_var%g'],...
    new_sigma^2*100);
if known_var 
    DCTO_locstr = strcat(DCTO_locstr,'_known') ; 
    KOH_locstr = strcat(KOH_locstr,'_known') ; 
    CTO_locstr = strcat(CTO_locstr,'_known') ; 
end
load(DCTO_locstr,'DCTO_results')
load(KOH_locstr,'KOH_results')
load(CTO_locstr,'CTO_results')
burn_in = DCTO_results.settings.burn_in;

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
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimum');
% title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
% set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
suptitle(sprintf('DCTO and KOH+CTO, obs. error var %g, discrep. %g',...
    new_sigma^2,discrep));
flushLegend(lg1,f1.Children(5),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    
%% Get posterior variance of point estimators of t1, t2
clc ; clearvars -except dpath ; close all ;

% Set number of times to perform calibration/design
n = 30;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, discrepancy version, whether modular
theta1 = 2;
% discrep = 0;
sigma = sqrt(0.05);
modular = false;
obs_discrep = true; % Whether or not to include discrep term for real obs
known_var = true; % Whether or not to treat new_sigma as known
if known_var, sigma2=sigma^2; else sigma2=0.05 ; end

% Set des_x size
des_x_size = 20;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Start loop for recording point estimates of t1,t2
for jj = 1:n
    
close all; disp(jj);

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2; mean_y = design.mean_y;
std_y = design.std_y;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'emulator',false,...
    'modular',modular,...
    'ObsVar',sigma2);

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
    'des_discrep',false,...
    'emulator',false,...
    'modular',modular,...
    'EmulatorMean',mean_sim,...
    'ObsVar',sigma2);

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
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim);

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
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
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
suptitle(sprintf('DCTO and KOH+CTO using observation error var %g',...
    sigma^2));
flushLegend(lg1,f1.Children(5),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');

% Record this loop's point estimates of t1,t2
DCTO_theta1s(jj) = mean(DCTO_results.theta1(burn_in:end,:));
DCTO_theta2s(jj) = mean(DCTO_results.theta2(burn_in:end,:));
KOH_theta1s(jj) = mean(KOH_results.theta1(burn_in:end,:));
CTO_theta2s(jj) = mean(CTO_results.theta1(burn_in:end,:));

end


% Create structure containing the point estimates and their variances
estimates = struct('DCTO_theta1s',DCTO_theta1s,...
    'DCTO_theta2s',DCTO_theta2s,'KOH_theta1s',KOH_theta1s,...
    'CTO_theta2s',CTO_theta2s,'DCTO_theta1_var',var(DCTO_theta1s),...
    'DCTO_theta2_var',var(DCTO_theta2s),...
    'KOH_theta1_var',var(KOH_theta1s),...
    'CTO_theta2_var',var(CTO_theta2s));

% Save it
discrep_str = int2str(discrep);
savestr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-07_Posterior_point_estimates_with_obs_noise-discrep',...
    discrep_str]);
% save(savestr,'estimates');

%% Examine all results sampling posterior point estimate variance
clc ; clearvars -except dpath ; close all ;

titles = {'No discrepancy:'    ;
    'Discrepancy form g_1 (small):'              ;
    'Discrepancy form g_1 (large):'              ;
    'Discrepancy form g_2 (small):'   ;
    'Discrepancy form g_2 (large):'   ;
    'Discrepancy form g_3 (small):'     ;
    'Discrepancy form g_3 (large):'}   ;

for discrep = 0:6
    discrep_str = int2str(discrep);
    loadstr = sprintf(['C:\\Users\\carle\\Documents',...
        '\\MATLAB\\NSF DEMS\\Phase 1\\',...
        'dual_calib\\dual_calib_stored_data\\'...
        '2019-08-07_Posterior_point_estimates_with_obs_noise-discrep',...
        discrep_str]);
    load(loadstr); 
    fields = {'DCTO_theta1s','DCTO_theta2s','KOH_theta1s','CTO_theta2s'};
%     fields={'DCTO_theta1s','DCTO_theta2s','KOH_theta1s','CTO_theta2s',...
%         'DCTO_theta1_var','DCTO_theta2_var','KOH_theta1_var',...
%         'CTO_theta2_var',};
    fprintf(strcat(titles{discrep+1},'\n'));
    disp(rmfield(estimates,fields));
%     disp(range(estimates.DCTO_theta1s));
%     disp(range(estimates.DCTO_theta2s));
%     disp(range(estimates.KOH_theta1s));
%     disp(range(estimates.CTO_theta2s));
%     disp(mean(estimates.CTO_theta2s));
%     disp(sqrt(estimates.DCTO_theta1_var));
%     disp(sqrt(estimates.DCTO_theta2_var));
%     disp(sqrt(estimates.KOH_theta1_var));
%     disp(sqrt(estimates.CTO_theta2_var));
end

%% Gather results using fewer target outcomes
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 3;
discrep = 3;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, whether modular
theta1 = 2;
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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

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

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular);

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

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
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
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Gather DCTO and KOH+CTO results using 0 'obs' error for CTO
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 20;
discrep = 4;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, whether modular
theta1 = 2;
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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

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

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular);

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
    mean(theta1_samps),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

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
    'EmulatorMean',mean_sim,...
    'ObsVar',0);

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = KOH_results.settings.burn_in; 
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = CTO_results.settings.burn_in; 
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Prepare to compare DCTO and CTO after KOH, mid-run
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 4;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .25 ; burn_in=M*b;

% Set real theta1, whether modular
theta1 = 2;
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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

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

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular);

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
    mean(theta1_samps),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

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

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
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
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Compare DCTO and KOH+CTO using no target discrepancy
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 4;

% Set number of draws, burn_in for each mcmc:
M = 8e2; b = .25 ; burn_in=M*b;

% Set real theta1, whether modular, covariance hyperparams
theta1 = 2;
modular = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
DesVar = 0.05 ; % error/tolerance
obs_rho_beta_params = [2,.4];
obs_lambda_gam_params = [10,10];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [150,4];

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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

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

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',DesVar,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params);

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
    'modular',false,...
    'EmulatorMean',mean_sim,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params);

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
    mean(theta1_samps),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

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
    'obs_discrep',des_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim,...
    'obs_var',DesVar,...
    'obs_rho_beta_params',des_rho_beta_params,...
    'obs_lambda_gam_params',des_lambda_gam_params);

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = KOH_results.settings.burn_in; 
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = CTO_results.settings.burn_in; 
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')


%% Compare DCTO with KOH+CTO including obs discrep info in the CTO
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 6;

% Set number of draws, burn_in for each mcmc:
M = 2e3; b = .5 ; burn_in=M*b;

% Set real theta1, whether modular, covariance hyperparams
theta1 = 2;
modular = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
des_var = 0; % target error/tolerance
obs_var = 0.05; % observation error
obs_rho_beta_params = [2,.4];
obs_lambda_gam_params = [10,10];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [100,4];

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
obs_x = design.obs_x; obs_t2 = design.obs_t2;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise / sigma * sqrt(obs_var);

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Load raw data for emulator, just to get mean and std
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-09-10_dual_calib_raw_data_for_emulator']);
mean_y = mean(sim_y) ; std_y = std(sim_y) ; 

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(a,b,c) dual_calib_example_fn(...
            a,xmin,xrange,b,t1min,t1range,c,t2min,t2range,...
            mean_y,std_y);

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params);

% Perform dual calibration
DCTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform KOH calibration

% Get "real" observations
% Get same observations used in DCTO version:
KOH_obs_y = DCTO_results.settings.obs_y * std_y + mean_y;

% Our x is now bivariate, incorporating the former x plus now also what was
% t2. This is because t2 is treated here as a control variable, not
% something over which we are attempting to optimize.

% Now set desired observations; since we're doing KOH, there are none
KOH_des_x = zeros(0,2); KOH_des_y = [];

% Define objective function
mean_sim_KOH = @(a,b,c) dual_calib_example_fn(...
    a(:,1),xmin,xrange,...
    b,t1min,t1range,... 
    a(:,2),t2min,t2range,...
    mean_y,std_y);

% Get settings
settings = MCMC_dual_calib_settings(zeros(0,2),sim_t1,[],sim_y,...
    [obs_x,obs_t2],[],KOH_obs_y,KOH_des_x,KOH_des_y,...
    'min_x',[xmin,t2min],...
    'range_x',[xrange,t2range],...
    'min_t1',t1min,'range_t1',t1range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim_KOH,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params);

% Perform calibration
KOH_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform CTO after KOH

% Get theta1 dist from KOH:
theta1_samps = ...
    (KOH_results.theta1(KOH_results.settings.burn_in:end) - ...
    KOH_results.settings.min_t1) / KOH_results.settings.range_t1 ;

% Define objective function
mean_sim_CTO = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,...
    mean(theta1_samps),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

% Now set desired observations
CTO_obs_x_01 = linspace(0,1,des_x_size)';
CTO_obs_x = CTO_obs_x_01 * xrange + xmin;
CTO_obs_y = zeros(size(CTO_obs_x,1),1);
obs_t2 = [] ; % We set obs_t2 to be empty, since we only have a 
% KOH calibration here, no DCTO secondary calibration. Similarly we will
% set des_x,des_y to be empty
des_x = [] ; des_y = [];

% Get additional observation discrepancy information learned from KOH
mean_obs = @(x,t) zeros(size(x));
add_nug = @(X) X+eye(size(X))*1e-4;
obs_discrep_vals = KOH_results.settings.obs_y - ...
    mean_sim(KOH_results.settings.obs_x(:,1),...
        repmat(mean(theta1_samps),...
            size(KOH_results.settings.obs_x,1),1),...
            KOH_results.settings.obs_x(:,2));
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    mean_obs(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1))))*...
    (y - mean_obs(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1)))) * ...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
additional_discrep_mean = @(t2) updated_mean(...
    obs_discrep_vals,...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));
additional_discrep_cov = @(t2) updated_cov(...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    CTO_obs_x,obs_t2,CTO_obs_y,des_x,des_y,...
    'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',des_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim_CTO,...
    'obs_var',des_var,...
    'additional_discrep_mean',additional_discrep_mean,...
    'additional_discrep_cov',additional_discrep_cov,...
    'obs_rho_beta_params',des_rho_beta_params,...
    'obs_lambda_gam_params',des_lambda_gam_params);

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = KOH_results.settings.burn_in; 
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = CTO_results.settings.burn_in; 
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimum');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-22_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-22_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-22_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Get target outcomes near pareto front 
clc ; clearvars -except dpath ; close all ;

% Choose discrepancy version, true theta1 value, # of targets
discrep = 6 ; 
theta1 = 2 ;
des_x_size = 15 ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Find optimal theta2 value
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);

% Get target observation locations
des_x = linspace(xmin,xmin+xrange,des_x_size)';

% Get target outcomes
des_y = ...
    dual_calib_example_fn(des_x,0,1,theta1,0,1,theta2,0,1,0,1,discrep) *.9;

% Save
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-21_nearby_target_outcomes_discrep',int2str(discrep),...
    '_des_x_size',int2str(des_x_size)]);
% save(locstr,'des_y');

%% Compare DCTO with KOH+CTO - nearby targets, updated discrep in CTO
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 6;

% Set number of draws, burn_in for each mcmc:
M = 2e3; b = .5 ; burn_in=M*b;

% Set real theta1, whether modular, covar hyperparams
theta1 = 2;
modular = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
des_var = 0 ; % target error/tolerance
obs_var = 0.05 ; % observation error
obs_rho_beta_params = [2,.4];
obs_lambda_gam_params = [10,10];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [100,4];

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

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise / sigma * sqrt(obs_var);

% Now load desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-21_nearby_target_outcomes_discrep',int2str(discrep),...
    '_des_x_size',int2str(des_x_size)]);
load(locstr,'des_y');

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(a,b,c) dual_calib_example_fn(...
            a,xmin,xrange,b,t1min,t1range,c,t2min,t2range,...
            mean_y,std_y);

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params,...
    'des_var',des_var,...
    'obs_var',obs_var);

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
KOH_des_x = zeros(0,2); KOH_des_y = [];

% Define objective function
mean_sim_KOH = @(a,b,c) dual_calib_example_fn(...
    a(:,1),xmin,xrange,...
    b,t1min,t1range,... 
    a(:,2),t2min,t2range,...
    mean_y,std_y);

% Get settings
settings = MCMC_dual_calib_settings(zeros(0,2),sim_t1,[],sim_y,...
    [obs_x,obs_t2],[],obs_y,KOH_des_x,KOH_des_y,'min_x',[xmin,t2min],...
    'range_x',[xrange,t2range],...
    'min_t1',t1min,'range_t1',t1range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim_KOH,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'obs_var',obs_var);

% Perform calibration
KOH_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform CTO after KOH

% Get theta1 dist from KOH:
theta1_samps = ...
    (KOH_results.theta1(KOH_results.settings.burn_in:end) - ...
    KOH_results.settings.min_t1) / KOH_results.settings.range_t1 ;

% Define objective function
mean_sim_CTO = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,...
    mean(theta1_samps),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

% Now set desired observations
CTO_obs_x_01 = linspace(0,1,des_x_size)' ;
CTO_obs_x = CTO_obs_x_01 * xrange + xmin;
CTO_obs_y = des_y;
CTO_obs_t2 = [] ; % We set obs_t2 to be empty, since we only have a 
% KOH calibration here, no DCTO secondary calibration. Similarly we will
% set des_x,des_y to be empty
des_x = [] ; des_y = [];

% Get additional observation discrepancy information learned from KOH
mean_obs = @(x,t) zeros(size(x));
add_nug = @(X) X+eye(size(X))*1e-4;
obs_discrep_vals = KOH_results.settings.obs_y - ...
    mean_sim(KOH_results.settings.obs_x(:,1),...
        repmat(mean(theta1_samps),...
            size(KOH_results.settings.obs_x,1),1),...
            KOH_results.settings.obs_x(:,2));
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    mean_obs(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1))))*...
    (y - mean_obs(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1)))) * ...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
additional_discrep_mean = @(t2) updated_mean(...
    obs_discrep_vals,...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));
additional_discrep_cov = @(t2) updated_cov(...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    CTO_obs_x,CTO_obs_t2,CTO_obs_y,des_x,des_y,...
    'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',des_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim_CTO,...
    'obs_var',des_var,...
    'obs_rho_beta_params',des_rho_beta_params,...
    'obs_lambda_gam_params',des_lambda_gam_params,...
    'additional_discrep_mean',additional_discrep_mean,...
    'additional_discrep_cov',additional_discrep_cov);

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = KOH_results.settings.burn_in; 
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = CTO_results.settings.burn_in; 
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Perform DCTO and KOH+CTO using emulator
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 6;

% Set number of draws, burn_in for each mcmc:
M = 2e3; b = .5 ; burn_in=M*b;

% Set real theta1, whether modular, covariance hyperparams
theta1 = 2;
modular = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
des_var = 0 ; % target error/tolerance
obs_var = 0.05 ; % observation error
obs_rho_beta_params = [2,.4];
obs_lambda_gam_params = [10,10];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [150,4];

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Load raw data for emulator
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-09-10_dual_calib_raw_data_for_emulator']);

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise / sigma * sqrt(obs_var);

% Now load desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-21_nearby_target_outcomes_discrep',int2str(discrep),...
    '_des_x_size',int2str(des_x_size)]);
load(locstr,'des_y');
des_y = des_y * 0 ; disp('Using target outcomes constant 0');

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x,1),1) ; 

% Emulator mean
mean_sim = @(a,b,c) zeros(size(a));

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',true,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params);

% TEMP set cov hyperparams to specified values
% settings.rho_proposal = @(r,s) ones(size(r)) * (.85 + .05*numel(r));
% settings.rho_prop_log_mh_correction = @(r_s,r) 0;
% settings.des_lambda_init = .6 ; settings.obs_lambda_init = 85;
% settings.lambda_proposal = @(lam,s) lam;
% settings.lambda_prop_log_mh_correction = @(lam_s,lam) 0;

% Perform dual calibration
DCTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform KOH calibration

% Get "real" observations
% Get same observations used in DCTO version:
obs_y = DCTO_results.settings.obs_y * std(sim_y) + mean(sim_y);

% Our x is now bivariate, incorporating the former x plus now also what was
% t2. This is because t2 is treated here as a control variable, not
% something over which we are attempting to optimize.

% Now set desired observations; since we're doing KOH, there are none
KOH_des_x = zeros(0,2); KOH_des_y = [];

% Define objective function
mean_sim_KOH = @(a,b,c) zeros(size(a,1),1);

% Get settings
settings = MCMC_dual_calib_settings([sim_x sim_t2],sim_t1,[],sim_y,...
    [obs_x,obs_t2],[],obs_y,KOH_des_x,KOH_des_y,'min_x',[xmin,t2min],...
    'range_x',[xrange,t2range],...
    'min_t1',t1min,'range_t1',t1range,'min_t2',[],'range_t2',[],...
    'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',true,...
    'modular',false,...
    'EmulatorMean',mean_sim_KOH,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params);

% Perform calibration
KOH_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform CTO after KOH

% Get theta1 dist from KOH:
theta1_samps = ...
    (KOH_results.theta1(KOH_results.settings.burn_in:end) - ...
    KOH_results.settings.min_t1) / KOH_results.settings.range_t1 ;

% Define objective function
mean_sim_CTO = @(a,b,c) zeros(size(a,1),1);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

% Now set desired observations
CTO_obs_x_01 = linspace(0,1,des_x_size)' ;
CTO_obs_x = CTO_obs_x_01 * xrange + xmin;
CTO_obs_y = des_y;
CTO_obs_t2 = mean(theta1_samps) * ones(size(CTO_obs_x)) * t1range + t1min ;
% We set obs_t2 to be the point estimate from KOH calibration
% We will set des_x,des_y to be empty
CTO_des_x = [] ; CTO_des_y = [];

% Get additional observation discrepancy information learned from KOH
mean_obs = @(x,t) zeros(size(x,1),1);
add_nug = @(X) X+eye(size(X))*1e-4;
obs_discrep_vals = KOH_results.settings.obs_y - ...
    mean_sim(KOH_results.settings.obs_x(:,1),...
        repmat(mean(theta1_samps),...
            size(KOH_results.settings.obs_x,1),1),...
            KOH_results.settings.obs_x(:,2));
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    mean_obs(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1))))*...
    (y - mean_obs(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1)))) * ...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
additional_discrep_mean = @(t2) updated_mean(...
    obs_discrep_vals,...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));
additional_discrep_cov = @(t2) updated_cov(...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));

% Get settings
settings = MCMC_dual_calib_settings([sim_x,sim_t1],sim_t2,[],sim_y,...
    [CTO_obs_x,CTO_obs_t2],[],CTO_obs_y,CTO_des_x,CTO_des_y,...
    'min_x',[xmin t1min],'range_x',[xrange t1range],...
    'min_t1',t2min,'range_t1',t2range,'min_t2',[],'range_t2',[],...
    'M',M,'burn_in',b,...
    'obs_discrep',des_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',true,...
    'modular',false,...
    'EmulatorMean',mean_sim_CTO,...
    'obs_var',des_var,...
    'additional_discrep_mean',additional_discrep_mean,...
    'additional_discrep_cov',additional_discrep_cov,...
    'obs_rho_beta_params',des_rho_beta_params,...
    'obs_lambda_gam_params',des_lambda_gam_params);

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = KOH_results.settings.burn_in; 
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = CTO_results.settings.burn_in; 
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.15);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.15);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Fiddle with covariance hyperparameter priors and observe the result
clc ; clearvars -except dpath ; close all ; 

% Set priors for covariance hyperparameters
lolp = @(lo) log( gampdf( lo, 150, 4 ) ) ;
% lolp = @(ld) log( gampdf( ld, 10,   10  ) ) ;
lorp = @(r ) sum( log( betapdf( r, 2, 0.4 ) ) ) ;
% lorp = @(r ) sum( log( betapdf( r, 2, 0.4 ) ) ) ;
ldlp = @(ld) log( gampdf( ld, 150,   4   ) ) ;
% ldlp = @(ld) log( gampdf( ld, 10,   10  ) ) ;
ldrp = @(r ) sum( log( betapdf( r, 2, 0.4 ) ) ) ;
% ldrp = @(r ) sum( log( betapdf( r, 2, 0.4 ) ) ) ;

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 6;

% Load nonmodular DCTO results
locstr = ['C:\\Users\\carle\\Documents\\MATLAB\\NSF DEMS\\'...
    'Phase 1\\dual_calib\\dual_calib_stored_data\\'...
    '2019-09-05_nonmodular_DCTO_discrep',int2str(discrep)];
load(locstr);

% Set number of draws, burn_in for each mcmc:
M = 2e3; b = .5 ; burn_in=M*b;

% Set real theta1, whether modular
theta1 = 2;
modular = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
des_var = 0 ; % target error/tolerance
obs_var = 0.05 ; % observation error

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Load raw data for emulator
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-08-28_dual_calib_raw_data_discrep',int2str(discrep)]);

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise / sigma * sqrt(obs_var);

% Now load desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-21_nearby_target_outcomes_discrep',int2str(discrep),...
    '_des_x_size',int2str(des_x_size)]);
load(locstr,'des_y');
des_y = des_y * 0 ; disp('Using target outcomes constant 0');

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x,1),1) ; 

% Emulator mean
mean_sim = @(a,b,c) zeros(size(a));

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',true,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var);

% set cov priors to specified values
settings.log_obs_lambda_prior = lolp ;
settings.log_obs_rho_prior = lorp ; 
settings.log_des_lambda_prior = ldlp ;
settings.log_des_rho_prior = ldrp ; 
settings.obs_lambda_init = 60;

% Perform dual calibration
DCTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = DCTO_results.settings.burn_in; 
histogram(nonmodres.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','Nonmod','Mod','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = nonmodres.settings.burn_in; 
histogram(nonmodres.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.15);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.15);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','Nonmod','Mod','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%%%%%% 
% Another figure: compare covariance hyperparameters
figure();
subplot(2,3,1);
histogram(DCTO_results.obs_rho(burn_in:end,1),'BinWidth',0.05,...
    'Normalization','pdf') ; 
hold on;
histogram(nonmodres.obs_rho(burn_in:end,1),'BinWidth',0.05,...
    'Normalization','pdf') ; 
title('\rho_{obs}_1');

subplot(2,3,2);
histogram(DCTO_results.obs_rho(burn_in:end,2),'BinWidth',0.05,...
    'Normalization','pdf') ; 
hold on;
histogram(nonmodres.obs_rho(burn_in:end,2),'BinWidth',0.05,...
    'Normalization','pdf') ; 
title('\rho_{obs}_2');

subplot(2,3,3);
histogram(DCTO_results.obs_lambda(burn_in:end),'BinWidth',10,...
    'Normalization','pdf') ; 
hold on;
histogram(nonmodres.obs_lambda(burn_in:end),'BinWidth',10,...
    'Normalization','pdf') ; 
title('\lambda_{obs}');

subplot(2,3,4);
histogram(DCTO_results.des_rho(burn_in:end),'BinWidth',0.02,...
    'Normalization','pdf') ; 
hold on;
histogram(nonmodres.des_rho(burn_in:end),'BinWidth',0.02,...
    'Normalization','pdf') ; 
title('\rho_{des}');

subplot(2,3,5);
histogram(DCTO_results.des_lambda(burn_in:end),'BinWidth',0.2,...
    'Normalization','pdf') ; 
hold on;
histogram(nonmodres.des_lambda(burn_in:end),'BinWidth',0.2,...
    'Normalization','pdf') ; 
title('\lambda_{des}');


%% Compare two versions of modularity
clc ; clearvars -except dpath ; close all ; 

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 5;

% Set number of draws, burn_in for each mcmc:
M = 2e3; b = .5;

% Set real theta1, whether emulator, whether modular, covar hyperparams
theta1 = 2;
emulator = false;
modular = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
des_var = 0 ; % target error/tolerance
obs_var = 0.05 ; % observation error
obs_rho_beta_params = [2,.4];
obs_lambda_gam_params = [10,10];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [150,4];

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Load raw data for emulator, get mean and std
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-09-10_dual_calib_raw_data_for_emulator']);
mean_y = mean(sim_y) ; std_y = std(sim_y);

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2;

% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn((obs_x-xmin)/xrange,...
    xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,...
    (obs_t2-t2min)/t2range,t2min,t2range,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise / sigma * sqrt(obs_var);

% Now load desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-21_nearby_target_outcomes_discrep',int2str(discrep),...
    '_des_x_size',int2str(des_x_size)]);
load(locstr,'des_y');
des_y = des_y * 0 ; disp('Using target outcomes constant 0');

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x,1),1) ; 

if emulator
    % Emulator mean
    mean_sim = @(a,b,c) zeros(size(a)); 
else
    % Since we are not using emulator, empty out simulator observations
    sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

    % Emulator mean
    mean_sim = @(a,b,c) dual_calib_example_fn(...
                a,xmin,xrange,b,t1min,t1range,c,t2min,t2range,...
                mean_y,std_y);
end

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'M',M,'burn_in',b,...
    'mean_y',mean_y,...
    'std_y',std_y,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',emulator,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params);

% Perform dual calibration
DCTO_fullmod_results = MCMC_dual_calib(settings);
DCTO_halfmod_results = MCMC_dual_calib_altmod(settings);


%%%%%%%%%%%%%%%%%
% Now make figures 

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
burn_in = DCTO_fullmod_results.settings.burn_in; 
histogram(DCTO_halfmod_results.theta1(burn_in+1:end),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
histogram(DCTO_fullmod_results.theta1(burn_in+1:end),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','Fullmod','Halfmod','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(DCTO_halfmod_results.theta2(burn_in+1:end,:),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.15);
histogram(DCTO_fullmod_results.theta2(burn_in+1:end,:),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.15);
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],t2min,t2min+t2range);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','Fullmod','Halfmod','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')


%% Explore objective function output when theta1 is some fn of theta2
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 4.5;
t1range = 9;
t2min = -3;
t2range = 2;

% Set discrepancy
discrep = 6;

% Set base theta1
low_theta1 = 1.5
high_theta1 = 2.25

% The function
% t1_fn = @(t2) t1min + t1range - t1range / t2range * (t2 - t2min) ;
% t1_fn = @(t2) t1min * 7/6 + (t1range - t1min/2) *...
%     exp(-2*(t2-(t2min+t2range*7/10)).^10) ;
t1_fn = @(x) low_theta1 + ...
    (high_theta1-low_theta1) * ...
    exp(40*(x-t2min)/t2range-20)./(1+exp(40*(x-t2min)/t2range-20));
% t1_fn = @(x) base_theta1 - ...
%     2.25 *exp(40*(x-t2min)/t2range-20)./(1+exp(40*(x-t2min)/t2range-20));
% t1_fn = @(x) high_theta1 - ...
%     (high_theta1 - low_theta1) * ...
%     exp(40*(x-t2min)/t2range-20)./(1+exp(40*(x-t2min)/t2range-20));
% t1_fn = @(x) low_theta1 + ...
%     (high_theta1-low_theta1) * ...
%     exp(20*((x-t2min)/t2range).^2-10)./...
%     (1+exp(20*((x-t2min)/t2range).^2-10));
% t1_fn = @(x) low_theta1 + ...
%     (high_theta1-low_theta1) * ...
%     exp(80*((x-t2min)/t2range)-40)./...
%     (1+exp(80*((x-t2min)/t2range)-40));
t1_fn = @(x) low_theta1 + ...
    (high_theta1-low_theta1) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));
t1_fn = @(x) high_theta1 - ...
    (high_theta1-low_theta1) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));
% t1_fn=@(x)2.5 * ones(size(x))

% Let's take a look at the objective function values for set x, using true
% t1 function as well as just base theta1
x = 1;
t2 = linspace(t2min,t2min+t2range,10000)';
t1 = t1_fn(t2);
y = dual_calib_example_fn(x,xmin,xrange,t1,t1min,t1range,...
    t2,t2min,t2range,0,1,discrep,false);
y_wrong = dual_calib_example_fn(x,xmin,xrange,low_theta1,t1min,t1range,...
    t2,t2min,t2range,0,1,discrep,false);
y_wrong2= dual_calib_example_fn(x,xmin,xrange,high_theta1,t1min,t1range,...
    t2,t2min,t2range,0,1,discrep,false);
figure('Position',[10 10 1100 300]);
subplot(1,3,1);
plot(t2,t1);
xlabel('t2');ylabel('t1');
subplot(1,3,2);
plot(t1,y,'.');
xlabel('t1');ylabel('y');
subplot(1,3,3);
plot(t2,y);
xlabel('t2');ylabel('y');
hold on;
plot(t2,y_wrong);
plot(t2,y_wrong2);
[m,i] = min(y) ; t2opt = t2(i)
t1opt = t1(i)

%% Get "observed" data when t1 is a function of t2
clc ; clearvars -except dpath ; close all ; 

% Set observation set size
obs_size = 30 ;

% Set discrepancy version
discrep = 6 ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% The function
% t1_fn = @(t2) t1min + t1range - t1range / t2range * (t2 - t2min) ;
t1_fn = @(t2) t1min * 7/6 + (t1range - t1min/2) *...
    exp(-50*(t2-(t2min+t2range*3/10)).^4) ;


% Get DoE
Obs = lhsdesign(obs_size,2) ; 
obs_x = Obs(:,1) * xrange + xmin ; obs_t2 = Obs(:,2) * t2range + t2min ; 
% Get t1 values
obs_t1 = t1_fn(obs_t2);

% Now get objective function values
obs_y = dual_calib_example_fn(obs_x,0,1,obs_t1,0,1,obs_t2,0,1,0,1,discrep);

% Now save it all
Obs = struct('obs_x',obs_x,'obs_t2',obs_t2,'obs_y_noiseless',obs_y);
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-09-13_state_aware_obs_discrep',...
    int2str(discrep)]);
% save(locstr,'Obs');

%% Compare DCTO and KOH+CTO in the state-aware case
clc ; clearvars -except dpath ; close all ;

% Set des_x size and discrepancy
des_x_size = 15;
discrep = 2;

% Set number of draws, burn_in for each mcmc:
M = 2e3; b = .5 ; burn_in=M*b;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Settings
modular = false;
informative_targets = true;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
des_var = 0; % target error/tolerance
obs_var = 0.05; % observation error
obs_rho_beta_params = [2,.4];
obs_lambda_gam_params = [10,10];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [100,4];

% True theta1 function
t1_fn = @(t2) t1min + t1range - t1range / t2range * (t2 - t2min) ;

% Load saved design
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-03-15_obs_design_and_sim_mean_std'],'design');
obs_x = design.obs_x; obs_t2 = design.obs_t2;

% Make a col vector based on true theta1
obs_t1 = t1_fn(obs_t2);

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn(obs_x,0,1,...
    obs_t1,0,1,...
    obs_t2,0,1,0,1,discrep);

% Now noise it up
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
% obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;
% Load noise from file, so all runs can use the same noise.
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\2019-07-22_obs_noise']);
obs_y = obs_y_noiseless + tempnoise / sigma * sqrt(obs_var);

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
if informative_targets    
    % Find true minimum at des_x
    fn=@(t2)dual_calib_example_fn(.5,0,1,t1_fn(t2),0,1,t2,0,1,0,1,discrep);
    theta2 = ...
        fmincon(fn,rand*t2range + t2min,[],[],[],[],t2min,t2min+t2range);
        true_min = dual_calib_example_fn(des_x,0,1,...
            t1_fn(theta2) * ones(size(des_x)),0,1,...
            theta2,0,1,0,1,discrep);
        des_y = true_min - 2*sqrt(obs_var) ; 
else
    des_y = zeros(size(des_x,1),1);
end
        

% Load raw data for emulator, just to get mean and std
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-09-10_dual_calib_raw_data_for_emulator']);
mean_y = mean(sim_y) ; std_y = std(sim_y) ; 

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
avg_disc=0 ; % Actually nevermind
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(a,b,c) dual_calib_example_fn(...
            a,xmin,xrange,b,t1min,t1range,c,t2min,t2range,...
            mean_y,std_y);

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params);

% Perform dual calibration
DCTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform KOH calibration

% Get "real" observations
% Get same observations used in DCTO version:
KOH_obs_y = DCTO_results.settings.obs_y * std_y + mean_y;

% Our x is now bivariate, incorporating the former x plus now also what was
% t2. This is because t2 is treated here as a control variable, not
% something over which we are attempting to optimize.

% Now set desired observations; since we're doing KOH, there are none
KOH_des_x = zeros(0,2); KOH_des_y = [];

% Define objective function
mean_sim_KOH = @(a,b,c) dual_calib_example_fn(...
    a(:,1),xmin,xrange,...
    b,t1min,t1range,... 
    a(:,2),t2min,t2range,...
    mean_y,std_y);

% Get settings
settings = MCMC_dual_calib_settings(zeros(0,2),sim_t1,[],sim_y,...
    [obs_x,obs_t2],[],KOH_obs_y,KOH_des_x,KOH_des_y,...
    'min_x',[xmin,t2min],...
    'range_x',[xrange,t2range],...
    'min_t1',t1min,'range_t1',t1range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim_KOH,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params);

% Perform calibration
KOH_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now perform CTO after KOH

% Get theta1 dist from KOH:
theta1_samps = ...
    (KOH_results.theta1(KOH_results.settings.burn_in:end) - ...
    KOH_results.settings.min_t1) / KOH_results.settings.range_t1 ;

% Define objective function
mean_sim_CTO = @(a,b,c) dual_calib_example_fn(...
    a,xmin,xrange,...
    mean(theta1_samps),t1min,t1range,... 
    b,t2min,t2range,...
    mean_y,std_y);
% mean(theta1_samps) 1/6 randsample(theta1_samps,1)

% Now set desired observations
CTO_obs_x_01 = linspace(0,1,des_x_size)';
CTO_obs_x = des_x;
CTO_obs_y = des_y;
obs_t2 = [] ; % We set obs_t2 to be empty, since we only have a 
% KOH calibration here, no DCTO secondary calibration. Similarly we will
% set des_x,des_y to be empty
CTO_des_x = [] ; CTO_des_y = [];

% Get additional observation discrepancy information learned from KOH
mean_obs = @(x,t) zeros(size(x));
add_nug = @(X) X+eye(size(X))*1e-4;
obs_discrep_vals = KOH_results.settings.obs_y - ...
    mean_sim(KOH_results.settings.obs_x(:,1),...
        repmat(mean(theta1_samps),...
            size(KOH_results.settings.obs_x,1),1),...
            KOH_results.settings.obs_x(:,2));
prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
    gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
    [xp ones(size(xp,1),1).*t2p],lambda,false);
updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
    mean_obs(xp,t2p) + ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1))))*...
    (y - mean_obs(xo,t2o)) ;
updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
    add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
    (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
        add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
            obs_var*eye(size(xo,1)))) * ...
    add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
additional_discrep_mean = @(t2) updated_mean(...
    obs_discrep_vals,...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));
additional_discrep_cov = @(t2) updated_cov(...
    KOH_results.settings.obs_x(:,1),...
    KOH_results.settings.obs_x(:,2),...
    CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
    mean(KOH_results.obs_rho(burn_in:end,:)),...
    mean(KOH_results.obs_lambda(burn_in:end,:)));

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    CTO_obs_x,obs_t2,CTO_obs_y,CTO_des_x,CTO_des_y,...
    'min_x',xmin,'range_x',xrange,...
    'min_t1',t2min,'range_t1',t2range,'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',des_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'des_discrep',false,...
    'emulator',false,...
    'modular',false,...
    'EmulatorMean',mean_sim_CTO,...
    'obs_var',des_var,...
    'additional_discrep_mean',additional_discrep_mean,...
    'additional_discrep_cov',additional_discrep_cov,...
    'obs_rho_beta_params',des_rho_beta_params,...
    'obs_lambda_gam_params',des_lambda_gam_params);

% Perform calibration
CTO_results = MCMC_dual_calib(settings);

%%%%%%%%%%%%%%%%%
% Now make figures 
% Get optimal theta2 and corresponding theta1 value
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,0,1,t1_fn(z),0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,t2min+t2range/2,[],[],[],[],t2min,t2min+t2range);
theta1 = t1_fn(theta2);
% Get prior and posterior theta1
f1 = figure('pos',[10 10 550 225]);
lcol = [218 165 32]/255 ; % Color for line
subplot(1,2,1);
% Plot prior
fill([t1min t1min + t1range t1min + t1range t1min],...
    [0 0 1/t1range 1/t1range],'g','EdgeColor','none');
xlim([t1min t1min + t1range]);
hold on;
% Get a histogram of theta1 with true value marked
burn_in = KOH_results.settings.burn_in; 
histogram(KOH_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta1(burn_in+1:end),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','KOH','DCTO','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
burn_in = CTO_results.settings.burn_in; 
histogram(CTO_results.theta1(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65);
burn_in = DCTO_results.settings.burn_in; 
histogram(DCTO_results.theta2(burn_in+1:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','CTO','DCTO','Optimum');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');
    

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-22_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
KOH_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-22_KOH_discrep',discrep_str,'_des_x_size',des_x_size_str]);
CTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-22_CTO_after_KOH_discrep',discrep_str,...
    '_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')
% save(KOH_locstr,'KOH_results')
% save(CTO_locstr,'CTO_results')

%% Perform DCTO in state-aware context with sequential DoE
clc ; clearvars -except dpath ; close all ;

% Set initial and final observation set sizes
obs_initial_size = 10;
obs_final_size = 20;

% Settings
modular = true;
informative_targets = false;
obs_discrep = true; % Whether or not to include discrep term for real obs
des_discrep = true;
obs_discrep_use_MLEs = false;
des_var = 0; % target error/tolerance
obs_var = 0.05; % observation error
des_x_size = 15;
discrep = 6;

% Set number of draws, burn_in for each mcmc:
M = 6e3; b = .5 ; burn_in=M*b;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Observation and target discrepancy covariance hyperparameters
obs_rho_beta_params = [8,1];
obs_lambda_gam_params = [8,4];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [40,40];

% True theta1 function
% t1_fn = @(t2) t1min + t1range - t1range / t2range * (t2 - t2min) ;
% t1_fn = @(t2) t1min * 7/6 + (t1range - t1min/2) *...
%     exp(-50*(t2-(t2min+t2range*3/10)).^4) ;
% t1_fn = @(t2) t1min * 7/6 + (t1range - t1min/2) *...
%     exp(-2*(t2-(t2min+t2range*3/10)).^100) ;
% t1_fn = @(t2) t1min * 7/6 + (t1range - t1min/2) *...
%     exp(-2*(t2-(t2min+t2range*7/10)).^10) ;
% t1_fn = @(x) 1.75 + ...
%     2.25 *exp(40*(x-t2min)/t2range-20)./(1+exp(40*(x-t2min)/t2range-20));
% t1_fn = @(x) 2 + ...
%     (1) * ...
%     exp(80*((x-t2min)/t2range).^2-40)./...
%     (1+exp(80*((x-t2min)/t2range).^2-40));
% t1_fn = @(x) 1.5 + ...
%     (2.25-1.5) * ...
%     exp(40*((x-t2min)/t2range)-20)./...
%     (1+exp(40*((x-t2min)/t2range)-20));
t1_fn = @(x) 2.25 - ...
    (2.25-1.5) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));

% Load raw data for emulator, if just to get mean and std
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-09-10_dual_calib_raw_data_for_emulator']);
mean_y = mean(sim_y) ; std_y = std(sim_y) ; 

% Get initial design
X = lhsdesign(obs_initial_size,2);
obs_x = X(:,1) * xrange + xmin ; 
obs_t2 = X(:,2) * t2range + t2min ;

% Make a col vector based on true theta1
obs_t1 = t1_fn(obs_t2);

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn(obs_x,xmin,xrange,...
    obs_t1,t1min,t1range,...
    obs_t2,t2min,t2range,0,1,discrep,false);

% Now noise it up (add N(0,obs_var) noise to STANDARDIZED observations)
obs_y = obs_y_noiseless+ randn(obs_initial_size,1) * sqrt(obs_var) * std_y;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
if informative_targets    
    % Find true minimum at des_x
    fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
        t1_fn(t2),t1min,t1range,t2,t2min,t2range,0,1,discrep,false);
    theta2 = ...
        fmincon(fn,rand*t2range + t2min,[],[],[],[],t2min,t2min+t2range);
        true_min = dual_calib_example_fn(des_x,xmin,xrange,...
            t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
            theta2,t2min,t2range,0,1,discrep,false);
        des_y = true_min - 2*sqrt(obs_var) ; 
else
    des_y = zeros(size(des_x,1),1);
end
        
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t)...
    dual_calib_example_fn(x,xmin,xrange,t1_fn(t),...
    t1min,t1range,t,t2min,t2range,mean_y,std_y,discrep,false)-...
    dual_calib_example_fn(x,xmin,xrange,t1_fn(t),t1min,t1range,...
    t,t2min,t2range,...
    mean_y,std_y,0,false);
avg_disc = integral2(int_fn,xmin,xmin+xrange,t2min,t2min+t2range) / ...
    (xrange * t2range) ; 
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(a,b,c) dual_calib_example_fn(...
            a,xmin,xrange,b,t1min,t1range,c,t2min,t2range,...
            mean_y,std_y,0,true);
        
% True phenomenon
true_phenomenon = @(a,c) dual_calib_example_fn(...
    a,xmin,xrange,(t1_fn(c*t2range+t2min)-t1min)./t1range,...
    t1min,t1range,c,t2min,t2range,...
    mean_y,std_y,discrep,true) + randn(size(a)) * sqrt(obs_var);

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params,...
    'obs_final_size',obs_final_size,...
    'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
    'true_phenomenon',true_phenomenon);

% Perform dual calibration
SDOE_results = MCMC_dual_calib(settings);

%%% Now do DCTO with all observations made ahead of time
% Get observations
X = lhsdesign(obs_final_size,2);
obs_x = X(:,1) * xrange + xmin ; 
obs_t2 = X(:,2) * t2range + t2min ;

% Make a col vector based on true theta1
obs_t1 = t1_fn(obs_t2);

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = dual_calib_example_fn(obs_x,xmin,xrange,...
    obs_t1,t1min,t1range,...
    obs_t2,t2min,t2range,0,1,discrep,false);

% Now noise it up (add N(0,obs_var) noise to STANDARDIZED observations)
obs_y = obs_y_noiseless+ randn(obs_final_size,1) * sqrt(obs_var) * std_y;

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params,...
    'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
    'obs_final_size',obs_final_size);

% Perform dual calibration
PDOE_results = MCMC_dual_calib(settings);

%%%%%%%%%%
% Take a look at the results

% Get optimal theta2 and corresponding theta1 value
% Get and plot true theta2
fmfn =@(z) dual_calib_example_fn(...
    .75,xmin,xrange,t1_fn(z),t1min,t1range,...
    z,t2min,t2range,0,1,discrep,false);
theta2 = fmincon(fmfn,t2min+t2range/2,[],[],[],[],t2min,t2min+.5*t2range);
theta1 = t1_fn(theta2);

% Get prior and posterior theta1
f1 = figure('pos',[10 10 550 225]);
lcol = [218 165 32]/255 ; % Color for line
subplot(1,2,1);
% Plot prior
fill([t1min t1min + t1range t1min + t1range t1min],...
    [0 0 1/t1range 1/t1range],'g','EdgeColor','none');
xlim([t1min t1min + t1range]);
hold on;
% Get a histogram of theta1 with true value marked
burn_in = SDOE_results.settings.burn_in; 
histogram(PDOE_results.theta1(burn_in+1:end),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
histogram(SDOE_results.theta1(burn_in+1:end),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','PDoE','SDoE','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(PDOE_results.theta2(burn_in+1:end,:),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.15);
histogram(SDOE_results.theta2(burn_in+1:end,:),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.15);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','PDoE','SDoE','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');

% Also take a look at the observation locations of SDoE
figure();
plot(SDOE_results.settings.obs_x(1:obs_initial_size),...
    SDOE_results.settings.obs_t2(1:obs_initial_size),'*');
hold on
plot(SDOE_results.settings.obs_x((obs_initial_size+1):end),...
    SDOE_results.settings.obs_t2((obs_initial_size+1):end),'*');
plot(PDOE_results.settings.obs_x,PDOE_results.settings.obs_t2,'*');

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')

%% Create new example function designed for state-aware case
clc ; clearvars -except dpath ; close all ;

% Set mins and ranges
t1min = 1.5; t1range = 3 ; 
t2min = 0; t2range = 5 ;

t1 = linspace(t1min,t1min+t1range);
t2 = linspace(t2min,t2min+t2range);

[T1,T2] = meshgrid(t1,t2);

exfn = @(t1,t2) 4 - 2.5 * exp(-((t1-2).^2+(t2-3.5).^2)) - ...
    3 * exp(-((t1-3).^2+(t2-1.5).^2)) ; 

y = exfn(T1(:),T2(:)) ; 
Y = reshape(y,length(t1),length(t2));

% surf(T1,T2,Y);
contour(T2,T1,Y);

%%% Now let's look at the results of taking this as a state-aware case
% Define the function linking t1 to t2
high_theta1 = 4.25 ; low_theta1 = 2 ; 
t1_fn = @(x) high_theta1 - ...
    (high_theta1-low_theta1) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));

% Put the function on the contour plot
hold on;
plot(t2,t1_fn(t2));

% Let's take a look at the objective function values for set x, using true
% t1 function as well as just base theta1
t2 = linspace(t2min,t2min+t2range)';
t1 = t1_fn(t2);
y = exfn(t1,t2);
y_low = exfn(low_theta1,t2);
y_high = exfn(high_theta1,t2);


figure('Position',[10 10 1100 300]);
subplot(1,3,1);
plot(t2,t1);
xlabel('t2');ylabel('t1');
subplot(1,3,2);
plot(t1,y,'.');
xlabel('t1');ylabel('y');
subplot(1,3,3);
plot(t2,y);
xlabel('t2');ylabel('y');
hold on;
plot(t2,y_low);
plot(t2,y_high);
[m,i] = min(y) ; t2opt = t2(i)
t1opt = t1(i)

%% Use new example fn designed for state-aware case in DCTO with SDoE
clc ; clearvars -except dpath ; close all ;

% Set initial and final observation set sizes
obs_initial_size = 10;
obs_final_size = 20;

% Settings
modular = true;
informative_targets = false;
obs_discrep = false; % Whether or not to include discrep term for real obs
des_discrep = true;
obs_discrep_use_MLEs = false;
des_var = 0; % target error/tolerance
obs_var = 0.05; % observation error
des_x_size = 15;
discrep = 0;

% Set number of draws, burn_in for each mcmc:
M = 6e3; b = .5 ;

% Define inputs mins and ranges 
xmin = .75;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Observation and target discrepancy covariance hyperparameters
obs_rho_beta_params = [8,1];
obs_lambda_gam_params = [8,4];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [40,40];

% True theta1 function
high_theta1 = 4.25 ; low_theta1 = 2 ; 
t1_fn = @(x) high_theta1 - ...
    (high_theta1-low_theta1) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));

% Define the computer model & truth (version giving nonstandardized output)
model_fn_ns = @(x,t1,t2) 4 - (x*xrange+xmin).* ...
    (...
    4 * exp(-(((t1*t1range+t1min)-2).^2+...
    ((t2*t2range+t2min)-3.5).^2)) - ...
    1 * exp(-(((t1*t1range+t1min)-3).^2+((t2*t2range+t2min)-1.5).^2))...
    ) ; 
true_phenomenon_ns = @(x,t2) model_fn_ns(x,...
    (t1_fn(t2*t2range+t2min)-t1min)./t1range,t2);

% Get mean and std using computer model, use to define standardized version
X=lhsdesign(1000,3); 
Y=model_fn_ns(X(:,1),X(:,2),X(:,3));
mean_y = mean(Y) ; std_y = std(Y) ;
model_fn = @(x,t1,t2) (model_fn_ns(x,t1,t2) - mean_y)./std_y;
true_phenomenon = @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

% Get initial design
X = lhsdesign(obs_initial_size,2);
obs_x = X(:,1) * xrange + xmin ; 
obs_t2 = X(:,2) * t2range + t2min ;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
    (obs_t2-t2min)./t2range);

% Now noise it up (add N(0,obs_var) noise to STANDARDIZED observations)
obs_y = obs_y_noiseless+ randn(obs_initial_size,1) * sqrt(obs_var) * std_y;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
if informative_targets    
    % Find true minimum at des_x
    fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
        t1_fn(t2),t1min,t1range,t2,t2min,t2range,0,1,discrep,false);
    theta2 = ...
        fmincon(fn,rand*t2range + t2min,[],[],[],[],t2min,t2min+t2range);
        true_min = dual_calib_example_fn(des_x,xmin,xrange,...
            t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
            theta2,t2min,t2range,0,1,discrep,false);
        des_y = true_min - 2*sqrt(obs_var) ; 
else
    des_y = -4 * ones(size(des_x,1),1);
end
        
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t) true_phenomenon(x,t) - ...
    model_fn(x,(t1_fn(t*t2range+t2min)-t1min)./t1range,t) ;
avg_disc = integral2(int_fn,0,1,0,1) ; 
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(a,b,c) model_fn(a,b,c) ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params,...
    'obs_final_size',obs_final_size,...
    'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
    'true_phenomenon',true_phenomenon);

% Perform dual calibration
SDOE_results = MCMC_dual_calib(settings);

%%% Now do DCTO with all observations made ahead of time
% Get observations
X = lhsdesign(obs_final_size,2);
obs_x = X(:,1) * xrange + xmin ; 
obs_t2 = X(:,2) * t2range + t2min ;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
    (obs_t2-t2min)./t2range);

% Now noise it up (add N(0,obs_var) noise to STANDARDIZED observations)
obs_y = obs_y_noiseless+ randn(obs_final_size,1) * sqrt(obs_var) * std_y;

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',xmin,'range_x',xrange,...
    'min_t1',t1min,'range_t1',t1range,'min_t2',t2min,'range_t2',t2range,...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',obs_discrep,...
    'obs_discrep_mean',obs_discrep_mean,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'modular',modular,...
    'des_discrep',des_discrep,...
    'des_var',des_var,...
    'obs_var',obs_var,...
    'obs_rho_beta_params',obs_rho_beta_params,...
    'obs_lambda_gam_params',obs_lambda_gam_params,...
    'des_rho_beta_params',des_rho_beta_params,...
    'des_lambda_gam_params',des_lambda_gam_params,...
    'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
    'obs_final_size',obs_final_size);

% Perform dual calibration
PDOE_results = MCMC_dual_calib(settings);

%%%%%%%%%%
% Take a look at the results

% Get optimal theta2 and corresponding theta1 value
% Get and plot true theta2
fmfn =@(z) true_phenomenon_ns(.5,(z-t2min)./t2range);
theta2 = fmincon(fmfn,t2min+t2range/2,[],[],[],[],t2min,t2min+t2range);
theta1 = t1_fn(theta2);

% Get prior and posterior theta1
f1 = figure('pos',[10 10 550 225]);
lcol = [218 165 32]/255 ; % Color for line
subplot(1,2,1);
% Plot prior
fill([t1min t1min + t1range t1min + t1range t1min],...
    [0 0 1/t1range 1/t1range],'g','EdgeColor','none');
xlim([t1min t1min + t1range]);
hold on;
% Get a histogram of theta1 with true value marked
burn_in = SDOE_results.settings.burn_in; 
histogram(PDOE_results.theta1(burn_in+1:end),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.075);
histogram(SDOE_results.theta1(burn_in+1:end),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.075);
% Plot true theta1
plot([theta1 theta1],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg1 = legend('Prior dist.','PDoE','SDoE','True value');
% title('Prior and posterior distributions of \theta_1');
xlabel('\theta_1');
set(f1,'color','white');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([t2min t2min + t2range t2min + t2range t2min],...
    [0 0 1/t2range 1/t2range],'g','EdgeColor','none');
xlim([t2min t2min + t2range]);
hold on;
% Get a histogram of theta2 with true value marked
histogram(PDOE_results.theta2(burn_in+1:end,:),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',.65,'BinWidth',0.15);
histogram(SDOE_results.theta2(burn_in+1:end,:),...
    'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.65,'BinWidth',0.15);
% yyaxis left ;
plot([theta2 theta2],get(gca,'YLim'),'--','Color',lcol,'LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','PDoE','SDoE','Optimal value');
xlabel('\theta_2');

% suptitle(['Prior and posterior distributions of ',...
%     '\theta_1 (left) and \theta_2 (right)']);
% suptitle('CTO setting \theta_1=2.25');
flushLegend(lg1,f1.Children(4),'northeast');
flushLegend(lg2,f1.Children(2),'northeast');

% Also take a look at the observation locations of SDoE
figure();
plot(SDOE_results.settings.obs_x(1:obs_initial_size),...
    SDOE_results.settings.obs_t2(1:obs_initial_size),'*');
hold on
plot(SDOE_results.settings.obs_x((obs_initial_size+1):end),...
    SDOE_results.settings.obs_t2((obs_initial_size+1):end),'*');
plot(PDOE_results.settings.obs_x,PDOE_results.settings.obs_t2,'*');

% Save results
discrep_str = int2str(discrep);
des_x_size_str = int2str(des_x_size);
DCTO_locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-08-09_DCTO_discrep',discrep_str,'_des_x_size',des_x_size_str]);
% save(DCTO_locstr,'DCTO_results')

%% Gather results for SDoE DCTO
clc ; clearvars -except dpath ; close all ; 

for discrep = 0:6 % For each version of the discrepancy

    % Turn on observation discrepancy in the model iff discrep ~=0
    if discrep == 0
        obs_discrep = false;
    else 
        obs_discrep = true;
    end
    
    % Number of chains
    m = 2;

    for ii = 1:m
        % Set initial and final observation set sizes
        obs_initial_size = 0;
        obs_final_size = 20;

        % Settings
        modular = true;
        informative_targets = false;
        des_discrep = false;
        obs_discrep_use_MLEs = false;
        des_var = 0.5 % target error/tolerance
        obs_var = 0.05; % observation error
        des_x_size = 15;
        doplot=true;

        % Set number of draws, burn_in for each mcmc:
        M = 6e3; b = .5 ; burn_in = M*b;

        % Define inputs mins and ranges 
        xmin = .5;
        xrange = .5;
        t1min = 1.5;
        t1range = 3;
        t2min = 0;
        t2range = 5;

        % Observation and target discrepancy covariance hyperparameters
        obs_rho_beta_params = [8,1];
        obs_lambda_gam_params = [8,4];
        des_rho_beta_params = [2,.4];
        des_lambda_gam_params = [40,40];

        % True theta1 function
        high_theta1 = 2.25 ; low_theta1 = 1.5 ; 
%         t1_fn = @(x) high_theta1 - ...
%             (high_theta1-low_theta1) * ...
%             exp(40*((x-t2min)/t2range)-20)./...
%             (1+exp(40*((x-t2min)/t2range)-20));
        t1_fn = @(x) high_theta1 - ...
            (high_theta1-low_theta1) * ...
            exp(40*((x-t2min)/t2range)-20)./...
            (1+exp(40*((x-t2min)/t2range)-20));


        % Define comp. model & truth (version w/ nonstandardized output)
        model_fn_ns = @(x,t1,t2) dual_calib_example_fn(...
            x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,...
            0,1,0,true) ; 
        true_phenomenon_ns = @(x,t2) dual_calib_example_fn(...
            x,xmin,xrange,...
            (t1_fn(t2*t2range+t2min)-t1min)./t1range,t1min,t1range,...
            t2,t2min,t2range,0,1,discrep,true);
            

        % Get mean and std using comp. model, define standardized version
        X=lhsdesign(1000,3); 
        Y=model_fn_ns(X(:,1),X(:,2),X(:,3));
        mean_y = mean(Y) ; std_y = std(Y) ;
        model_fn = @(x,t1,t2) (model_fn_ns(x,t1,t2) - mean_y)./std_y;
        true_phenomenon = ...
            @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

        % Get initial design
        X = lhsdesign(obs_initial_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless+ ...
            randn(obs_initial_size,1) * sqrt(obs_var) * std_y;

        % Now set desired observations
        des_x = linspace(0,1,des_x_size)' * xrange + xmin;
        if informative_targets    
            % Find true minimum at des_x
            fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
                t1_fn(t2),t1min,t1range,...
                t2,t2min,t2range,0,1,discrep,false);
            theta2 = ...
                fmincon(fn,rand*t2range + t2min,...
                [],[],[],[],t2min,t2min+t2range);
                true_min = dual_calib_example_fn(des_x,xmin,xrange,...
                    t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
                    theta2,t2min,t2range,0,1,discrep,false);
                des_y = true_min - 2*sqrt(obs_var) ; 
        else
            des_y = 0 * ones(size(des_x,1),1);
        end

        % Since we are not using emulator, empty out simulator observations
        sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

        % And get the average discrepancy value to use as the mean
        int_fn =@(x,t) true_phenomenon(x,t) - ...
            model_fn(x,(t1_fn(t*t2range+t2min)-t1min)./t1range,t) ;
        avg_disc = integral2(int_fn,0,1,0,1) ; 
        fprintf('Average observation discrepancy: %g\n',avg_disc);
        obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

        % Emulator mean
        mean_sim = @(a,b,c) model_fn(a,b,c) ; 

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,...
            des_x,des_y,'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_final_size',obs_final_size,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'true_phenomenon',true_phenomenon,...
            'doplot',doplot);

        % Perform dual calibration
        % We need a loop because sometimes an ill-conditioned matrix early
        % in the burn-in makes the whole calibration fail.
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        
        sdoe_results.theta1(:,ii) = res.theta1 ;
        sdoe_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        sdoe_results.theta2(:,ii) = res.theta2 ;
        sdoe_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        sdoe_results.obs_rho(:,:,ii) = res.obs_rho ;
        sdoe_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        sdoe_results.obs_lambda(:,ii) = res.obs_lambda ;
        sdoe_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        sdoe_results.des_rho(:,ii) = res.des_rho ;
        sdoe_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        sdoe_results.des_lambda(:,ii) = res.des_lambda ;
        sdoe_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        sdoe_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        sdoe_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        sdoe_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        sdoe_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        sdoe_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        sdoe_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        sdoe_results.settings{ii} = settings ; 
        
        %%% Now do DCTO with all observations made ahead of time
        % Get observations
        X = lhsdesign(obs_final_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Make a col vector based on true theta1
        obs_t1 = t1_fn(obs_t2);

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless + ...
            randn(obs_final_size,1) * sqrt(obs_var) * std_y;

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,des_x,des_y,...
            'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'obs_final_size',obs_final_size,...
            'doplot',doplot);

        % Perform dual calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        pdoe_results.theta1(:,ii) = res.theta1 ;
        pdoe_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        pdoe_results.theta2(:,ii) = res.theta2 ;
        pdoe_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        pdoe_results.obs_rho(:,:,ii) = res.obs_rho ;
        pdoe_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        pdoe_results.obs_lambda(:,ii) = res.obs_lambda ;
        pdoe_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        pdoe_results.des_rho(:,ii) = res.des_rho ;
        pdoe_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        pdoe_results.des_lambda(:,ii) = res.des_lambda ;
        pdoe_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        pdoe_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        pdoe_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        pdoe_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        pdoe_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        pdoe_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        pdoe_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        pdoe_results.settings{ii} = settings ; 
        
        both_results.sdoe_results = sdoe_results;
        both_results.pdoe_results = pdoe_results;
        
        save('temp','both_results');
        
        % Close windows every now and again
        if mod(ii,10) == 0, close all ; end
        
        % Update location in the loop
        fprintf('\n####################################\n');
        fprintf('COMPLETED STEP %g/%g OF DISCREP %g\n',ii,m,discrep);
        fprintf('####################################\n');
        
        
    end
    
    % Get optimal theta2 and corresponding theta1 value
    % Get and plot true theta2
    fmfn =@(z) dual_calib_example_fn(...
        .75,xmin,xrange,t1_fn(z),t1min,t1range,...
        z,t2min,t2range,0,1,discrep,false);
    % Try two different start locations
    [theta2_1, fval_1] = fmincon(fmfn,t2min+t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [theta2_2, fval_2] = fmincon(fmfn,t2min+3*t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [~,idx] = min([fval_1 fval_2]);
    theta2s = [theta2_1 theta2_2] ; theta2 = theta2s(idx) ; 
    theta1 = t1_fn(theta2);
    
    % Add true parameter values, discrepancy version to results
    sdoe_results.discrepancy = discrep;
    pdoe_results.discrepancy = discrep;
    sdoe_results.true_theta1 = theta1;
    pdoe_results.true_theta1 = theta1;
    sdoe_results.true_theta2 = theta2;
    pdoe_results.true_theta2 = theta2;
    
    % Add results to full set of results
    results{discrep+1,1} = sdoe_results;
    results{discrep+1,2} = pdoe_results;
    
    % Save
    locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-10-04_SDOE_results_nobs' int2str(obs_initial_size) '-'...
    int2str(obs_final_size)]);
%     save(locstr,'results');
    
    % Close all those windows
    close all ;
    
    % Update location in the loop
    fprintf('####################################\n');
    fprintf('\nCOMPLETED DISCREPANCY VERSION %g\n',discrep);
    fprintf('####################################\n');
    
end


%% Collect results for DCTO vs KOH+CTO
clc ; clearvars -except dpath ; close all ; 

for discrep = 0:6 % For each version of the discrepancy

    % Turn on observation discrepancy in the model iff discrep ~=0
    if discrep == 0
        obs_discrep = false;
    else 
        obs_discrep = true;
    end
    
    % Number of chains
    m = 2;

    for ii = 1:m
        % Set observation set size
        obs_size = 20;

        % Settings
        modular = true;
        informative_targets = false;
        des_discrep = false;
        obs_discrep_use_MLEs = false;
        des_var = 0.5 % target error/tolerance
        obs_var = 0.05; % observation error
        des_x_size = 15;
        doplot = true;

        % Set number of draws, burn_in for each mcmc:
        M = 4e3; b = .5 ; burn_in = M*b;

        % Define inputs mins and ranges 
        xmin = .5;
        xrange = .5;
        t1min = 1.5;
        t1range = 3;
        t2min = 0;
        t2range = 5;

        % Observation and target discrepancy covariance hyperparameters
        obs_rho_beta_params = [8,1];
        obs_lambda_gam_params = [8,4];
        des_rho_beta_params = [2,.4];
        des_lambda_gam_params = [40,40];

        % True theta1 
        theta1 = 2 ;

        % Define comp. model & truth (version w/ nonstandardized output)
        model_fn_ns = @(x,t1,t2) dual_calib_example_fn(...
            x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,...
            0,1,0,true) ; 
        true_phenomenon_ns = @(x,t2) dual_calib_example_fn(...
            x,xmin,xrange,...
            (theta1-t1min)./t1range,t1min,t1range,...
            t2,t2min,t2range,0,1,discrep,true);
            

        % Get mean and std using comp. model, define standardized version
        X=lhsdesign(1000,3); 
        Y=model_fn_ns(X(:,1),X(:,2),X(:,3));
        mean_y = mean(Y) ; std_y = std(Y) ;
        model_fn = @(x,t1,t2) (model_fn_ns(x,t1,t2) - mean_y)./std_y;
        true_phenomenon = ...
            @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

        % Get initial design
        X = lhsdesign(obs_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless+ ...
            randn(obs_size,1) * sqrt(obs_var) * std_y;

        % Now set desired observations
        des_x = linspace(0,1,des_x_size)' * xrange + xmin;
        if informative_targets    
            % Find true minimum at des_x
            fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
                t1_fn(t2),t1min,t1range,...
                t2,t2min,t2range,0,1,discrep,false);
            theta2 = ...
                fmincon(fn,rand*t2range + t2min,...
                [],[],[],[],t2min,t2min+t2range);
                true_min = dual_calib_example_fn(des_x,xmin,xrange,...
                    t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
                    theta2,t2min,t2range,0,1,discrep,false);
                des_y = true_min - 2*sqrt(obs_var) ; 
        else
            des_y = 0 * ones(size(des_x,1),1);
        end

        % Since we are not using emulator, empty out simulator observations
        sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

        % And get the average discrepancy value to use as the mean
        int_fn =@(x,t) true_phenomenon(x,t) - ...
            model_fn(x,(theta1-t1min)./t1range,t) ;
%         avg_disc = integral2(int_fn,0,1,0,1) ; 
        avg_disc = mean(...
            true_phenomenon((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range) - ...
            model_fn((obs_x-xmin)./xrange,(theta1-t1min)./t1range,...
            (obs_t2-t2min)./t2range));
        fprintf('Average observation discrepancy: %g\n',avg_disc);
        obs_discrep_mean = @(x,t) avg_disc * ones(size(x,1),1) ; 
%         % Do this by linear regression on x,t2, without interaction.
%         obs_discrep_coeffs = ...
%             [ones(size(obs_x)) (obs_x-xmin)./xrange ...
%             (obs_t2-t2min)./t2range]\(obs_y-mean_y)./std_y ; 
%         obs_discrep_mean = @(x,t) ...
%             [ones(size(x,1),1) x t] * obs_discrep_coeffs ; 

        % Emulator mean
        mean_sim = @(a,b,c) model_fn(a,b,c) ; 

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,...
            des_x,des_y,'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'doplot',doplot,...
            'verbose',false);%,...
%             'des_var_est',true,...
%             'obs_var_est',true);

        % Perform dual calibration
        % We need a loop because sometimes an ill-conditioned matrix early
        % in the burn-in makes the whole calibration fail.
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        
        dcto_results.theta1(:,ii) = res.theta1 ;
        dcto_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        dcto_results.theta2(:,ii) = res.theta2 ;
        dcto_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        dcto_results.obs_rho(:,:,ii) = res.obs_rho ;
        dcto_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        dcto_results.obs_lambda(:,ii) = res.obs_lambda ;
        dcto_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        dcto_results.des_rho(:,ii) = res.des_rho ;
        dcto_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        dcto_results.des_lambda(:,ii) = res.des_lambda ;
        dcto_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        dcto_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        dcto_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        dcto_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        dcto_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        dcto_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        dcto_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        dcto_results.settings{ii} = settings ; 
        
        %%% Now do KOH calibration
        
        % Set observations
        KOH_obs_x = [obs_x obs_t2]; KOH_obs_t2 = [];
        
        % Set target outcomes to be empty
        KOH_des_x = zeros(0,2); KOH_des_y = [] ; 
        
        % Set sim inputs
        KOH_sim_x = zeros(0,2) ; % Avoids index error
        KOH_sim_t1 = [] ; KOH_sim_t2 = [] ; KOH_sim_y = [];
        
        % Set simulator mean for KOH
        KOH_mean_sim = @(x_t2,t1,junk) model_fn(x_t2(:,1),t1,x_t2(:,2));
        
        % Set obs discrep function for KOH
        KOH_obs_discrep_mean = @(x,t) ...
            obs_discrep_mean(x(:,1),x(:,2)) ;

        % Get settings for KOH
        settings = MCMC_dual_calib_settings(KOH_sim_x,KOH_sim_t1,...
            KOH_sim_t2,KOH_sim_y,...
            KOH_obs_x,KOH_obs_t2,obs_y,KOH_des_x,KOH_des_y,...
            'min_x',[xmin t2min],'range_x',[xrange t2range],...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',KOH_obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',KOH_mean_sim,...
            'modular',false,...
            'des_discrep',false,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'doplot',doplot,...
            'verbose',false);

        % Perform KOH calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        % Extract results
        khtc_results.theta1(:,ii) = res.theta1 ;
        khtc_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        khtc_results.obs_rho(:,:,ii) = res.obs_rho ;
        khtc_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        khtc_results.obs_lambda(:,ii) = res.obs_lambda ;
        khtc_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        khtc_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        khtc_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        khtc_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        khtc_results.settings{ii} = settings ; 
        
        %%%%%%
        % Now do CTO design using KOH calibration results
        
        % Set CTO obs to be target outcomes
        CTO_obs_x = des_x ; CTO_obs_t2 = [] ; CTO_obs_y = des_y; 
        CTO_obs_x_01 = (des_x - xmin)./xrange ; % Will come in handy
        
        % Set CTO targets to be empty
        CTO_des_x = [] ; CTO_des_y = [];
        
        % Set CTO "obs" discrepancy (actually target discrepancy)
        CTO_obs_discrep_mean = @(x,t) 0 * ones(size(x)) ; 
            
        % Get theta1 dist from KOH:
        theta1_samps =(khtc_results.theta1(burn_in:end,ii)-t1min)./t1range;

        % Define objective function
        CTO_mean_sim = @(x,t2,junk) model_fn(x, mean(theta1_samps), t2) ;

        % Get additional observation discrepancy informatn learned from KOH
        mean_obs = obs_discrep_mean;
        add_nug = @(X) X+eye(size(X))*1e-4;
        obs_discrep_vals = res.settings.obs_y - ...
            mean_sim(res.settings.obs_x(:,1),...
                repmat(mean(theta1_samps),...
                    size(res.settings.obs_x,1),1),...
                    res.settings.obs_x(:,2));
        prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
            gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
            [xp ones(size(xp,1),1).*t2p],lambda,false);
        updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda) ...
            mean_obs(xp,t2p) + ...
            (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
                add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
                    obs_var*eye(size(xo,1))))*...
            (y - mean_obs(xo,t2o)) ;
        updated_cov = @(xo,t2o,xp,t2p,rho,lambda) ...
            add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
            (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
                add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
                    obs_var*eye(size(xo,1)))) * ...
            add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
        additional_discrep_mean = @(t2) updated_mean(...
            obs_discrep_vals,...
            res.settings.obs_x(:,1),...
            res.settings.obs_x(:,2),...
            CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
            mean(res.obs_rho(burn_in:end,:)),...
            mean(res.obs_lambda(burn_in:end,:)));
        additional_discrep_cov = @(x,t2) updated_cov(...
            res.settings.obs_x(:,1),...
            res.settings.obs_x(:,2),...
            CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
            mean(res.obs_rho(burn_in:end,:)),...
            mean(res.obs_lambda(burn_in:end,:)));

        % Get settings
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            CTO_obs_x,CTO_obs_t2,CTO_obs_y,CTO_des_x,CTO_des_y,...
            'min_x',xmin,'range_x',xrange,...
            'min_t1',t2min,'range_t1',t2range,'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',des_discrep,...
            'obs_discrep_mean',CTO_obs_discrep_mean,...
            'des_discrep',false,...
            'emulator_use',false,...
            'modular',false,...
            'EmulatorMean',CTO_mean_sim,...
            'obs_var',des_var,...
            'additional_discrep_mean',additional_discrep_mean,...
            'additional_discrep_cov',additional_discrep_cov,...
            'obs_rho_beta_params',des_rho_beta_params,...
            'obs_lambda_gam_params',des_lambda_gam_params,...
            'doplot',doplot,...
            'verbose',false);

        % Perform CTO calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        % Extract results
        ctod_results.theta2(:,ii) = res.theta1 ;
        ctod_results.theta2_hat(ii) = mean(res.theta1(burn_in:end));
        ctod_results.des_rho(:,ii) = res.obs_rho ;
        ctod_results.des_rho_hat(ii) = mean(res.obs_rho(burn_in:end));
        ctod_results.des_lambda(:,ii) = res.obs_lambda ;
        ctod_results.des_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));        
        ctod_results.theta2_var(ii) = var(res.theta1(burn_in:end)) ;
        ctod_results.des_rho_var(ii) = ...
            var(res.obs_rho(burn_in:end)) ;
        ctod_results.des_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        ctod_results.settings = settings;

        
        all_results.dcto_results = dcto_results;
        all_results.khtc_results = khtc_results;
        all_results.ctod_results = ctod_results;
        
        save('temp','all_results');
        
        % Close windows every now and again
        if mod(ii,5) == 0, close all ; end
        
        % Update location in the loop
        fprintf('####################################\n');
        fprintf('COMPLETED STEP %g/%g OF DISCREP %g\n',ii,m,discrep);
        fprintf('####################################\n');
        
    end
    
    % Get optimal theta2 and corresponding theta1 value
    % Get and plot true theta2
    fmfn =@(z) dual_calib_example_fn(...
        .75,xmin,xrange,theta1,t1min,t1range,...
        z,t2min,t2range,0,1,discrep,false);
    % Try two different start locations
    [theta2_1, fval_1] = fmincon(fmfn,t2min+t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [theta2_2, fval_2] = fmincon(fmfn,t2min+3*t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [~,idx] = min([fval_1 fval_2]);
    theta2s = [theta2_1 theta2_2] ; theta2 = theta2s(idx) ; 
    
    % Add true parameter values, discrepancy version to results
    dcto_results.discrepancy = discrep;
    khtc_results.discrepancy = discrep;
    ctod_results.discrepancy = discrep;
    dcto_results.true_theta1 = theta1;
    khtc_results.true_theta1 = theta1;
    ctod_results.true_theta1 = theta1;
    dcto_results.true_theta2 = theta2;
    khtc_results.true_theta2 = theta2;
    ctod_results.true_theta2 = theta2;
    
    % Add results to full set of results
    results{discrep+1,1} = dcto_results;
    results{discrep+1,2} = khtc_results;
    results{discrep+1,3} = ctod_results;
    
    % Save
    locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-10-04_DCTO_vs_KOHCTO_results_temp']);
%     save(locstr,'results');
    
    % Close all those windows
    close all ;
    
    % Update location in the loop
    fprintf('####################################\n');
    fprintf('COMPLETED DISCREPANCY VERSION %g\n',discrep);
    fprintf('####################################\n');
    
end

%% Collect results for DCTO vs KOH+CTO w/ estimated target err var
clc ; clearvars -except dpath ; close all ; 

for discrep = 0:6 % For each version of the discrepancy

    % Turn on observation discrepancy in the model iff discrep ~=0
    if discrep == 0
        obs_discrep = false;
    else 
        obs_discrep = true;
    end
    
    % Number of chains
    m = 30;

    for ii = 1:m
        % Set observation set size
        obs_size = 20;

        % Settings
        modular = true;
        informative_targets = false;
        des_discrep = false;
        obs_discrep_use_MLEs = false;
        des_var = 0.5; % target error/tolerance
        obs_var = 0.05; % observation error
        des_x_size = 15;
        doplot = false;
        verbose = false;

        % Set number of draws, burn_in for each mcmc:
        M = 6e3; b = .5 ; burn_in = M*b;

        % Define inputs mins and ranges 
        xmin = .5;
        xrange = .5;
        t1min = 1.5;
        t1range = 3;
        t2min = 0;
        t2range = 5;

        % Observation and target discrepancy covariance hyperparameters
        obs_rho_beta_params = [8,1];
        obs_lambda_gam_params = [8,4];
        des_rho_beta_params = [2,.4];
        des_lambda_gam_params = [40,40];

        % True theta1 
        theta1 = 2 ;

        % Define comp. model & truth (version w/ nonstandardized output)
        model_fn_ns = @(x,t1,t2) dual_calib_example_fn(...
            x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,...
            0,1,0,true) ; 
        true_phenomenon_ns = @(x,t2) dual_calib_example_fn(...
            x,xmin,xrange,...
            (theta1-t1min)./t1range,t1min,t1range,...
            t2,t2min,t2range,0,1,discrep,true);
            

        % Get mean and std using comp. model, define standardized version
        X=lhsdesign(1000,3); 
        Y=model_fn_ns(X(:,1),X(:,2),X(:,3));
        mean_y = mean(Y) ; std_y = std(Y) ;
        model_fn = @(x,t1,t2) (model_fn_ns(x,t1,t2) - mean_y)./std_y;
        true_phenomenon = ...
            @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

        % Get initial design
        X = lhsdesign(obs_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless+ ...
            randn(obs_size,1) * sqrt(obs_var) * std_y;

        % Now set desired observations
        des_x = linspace(0,1,des_x_size)' * xrange + xmin;
        if informative_targets    
            % Find true minimum at des_x
            fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
                t1_fn(t2),t1min,t1range,...
                t2,t2min,t2range,0,1,discrep,false);
            theta2 = ...
                fmincon(fn,rand*t2range + t2min,...
                [],[],[],[],t2min,t2min+t2range);
                true_min = dual_calib_example_fn(des_x,xmin,xrange,...
                    t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
                    theta2,t2min,t2range,0,1,discrep,false);
                des_y = true_min - 2*sqrt(obs_var) ; 
        else
            des_y = 0 * ones(size(des_x,1),1);
        end

        % Since we are not using emulator, empty out simulator observations
        sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

        % And get the average discrepancy value to use as the mean
        int_fn =@(x,t) true_phenomenon(x,t) - ...
            model_fn(x,(theta1-t1min)./t1range,t) ;
%         avg_disc = integral2(int_fn,0,1,0,1) ; 
        avg_disc = mean(...
            true_phenomenon((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range) - ...
            model_fn((obs_x-xmin)./xrange,(theta1-t1min)./t1range,...
            (obs_t2-t2min)./t2range));
        fprintf('Average observation discrepancy: %g\n',avg_disc);
        obs_discrep_mean = @(x,t) avg_disc * ones(size(x,1),1) ; 
%         % Do this by linear regression on x,t2, without interaction.
%         obs_discrep_coeffs = ...
%             [ones(size(obs_x)) (obs_x-xmin)./xrange ...
%             (obs_t2-t2min)./t2range]\(obs_y-mean_y)./std_y ; 
%         obs_discrep_mean = @(x,t) ...
%             [ones(size(x,1),1) x t] * obs_discrep_coeffs ; 

        % Emulator mean
        mean_sim = @(a,b,c) model_fn(a,b,c) ; 

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,...
            des_x,des_y,'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'doplot',doplot,...
            'verbose',verbose,...
            'des_var_est',true,...
            'obs_var_est',false);

        % Perform dual calibration
        % We need a loop because sometimes an ill-conditioned matrix early
        % in the burn-in makes the whole calibration fail.
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        
        dcto_results.theta1(:,ii) = res.theta1 ;
        dcto_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        dcto_results.theta2(:,ii) = res.theta2 ;
        dcto_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        dcto_results.obs_rho(:,:,ii) = res.obs_rho ;
        dcto_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        dcto_results.obs_lambda(:,ii) = res.obs_lambda ;
        dcto_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        dcto_results.des_rho(:,ii) = res.des_rho ;
        dcto_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        dcto_results.des_lambda(:,ii) = res.des_lambda ;
        dcto_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
%         dcto_results.obs_var(:,ii) = res.obs_var;
        dcto_results.des_var(:,ii) = res.des_var;
%         dcto_results.obs_var_hat(ii) = mean(res.obs_var(burn_in:end));
        dcto_results.des_var_hat(ii) = mean(res.des_var(burn_in:end));
        dcto_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        dcto_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        dcto_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        dcto_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        dcto_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        dcto_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        dcto_results.settings{ii} = settings ; 
        
        %%% Now do KOH calibration
        
        % Set observations
        KOH_obs_x = [obs_x obs_t2]; KOH_obs_t2 = [];
        
        % Set target outcomes to be empty
        KOH_des_x = zeros(0,2); KOH_des_y = [] ; 
        
        % Set sim inputs
        KOH_sim_x = zeros(0,2) ; % Avoids index error
        KOH_sim_t1 = [] ; KOH_sim_t2 = [] ; KOH_sim_y = [];
        
        % Set simulator mean for KOH
        KOH_mean_sim = @(x_t2,t1,junk) model_fn(x_t2(:,1),t1,x_t2(:,2));
        
        % Set obs discrep function for KOH
        KOH_obs_discrep_mean = @(x,t) ...
            obs_discrep_mean(x(:,1),x(:,2)) ;

        % Get settings for KOH
        settings = MCMC_dual_calib_settings(KOH_sim_x,KOH_sim_t1,...
            KOH_sim_t2,KOH_sim_y,...
            KOH_obs_x,KOH_obs_t2,obs_y,KOH_des_x,KOH_des_y,...
            'min_x',[xmin t2min],'range_x',[xrange t2range],...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',KOH_obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',KOH_mean_sim,...
            'modular',false,...
            'des_discrep',false,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',false);

        % Perform KOH calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        % Extract results
        khtc_results.theta1(:,ii) = res.theta1 ;
        khtc_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        khtc_results.obs_rho(:,:,ii) = res.obs_rho ;
        khtc_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        khtc_results.obs_lambda(:,ii) = res.obs_lambda ;
        khtc_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        khtc_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        khtc_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        khtc_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
%         khtc_results.obs_var(:,ii) = res.obs_var;
%         khtc_results.obs_var_hat(ii) = mean(res.obs_var(burn_in:end));
        khtc_results.settings{ii} = settings ; 
        
        %%%%%%
        % Now do CTO design using KOH calibration results
        
        % Set CTO obs to be target outcomes
        CTO_obs_x = des_x ; CTO_obs_t2 = [] ; CTO_obs_y = des_y; 
        CTO_obs_x_01 = (des_x - xmin)./xrange ; % Will come in handy
        
        % Set CTO targets to be empty
        CTO_des_x = [] ; CTO_des_y = [];
        
        % Set CTO "obs" discrepancy (actually target discrepancy)
        CTO_obs_discrep_mean = @(x,t) 0 * ones(size(x)) ; 
            
        % Get theta1 dist from KOH:
        theta1_samps =(khtc_results.theta1(burn_in:end,ii)-t1min)./t1range;

        % Define objective function
        CTO_mean_sim = @(x,t2,junk) model_fn(x, mean(theta1_samps), t2) ;

        % Get additional observation discrepancy informatn learned from KOH
        mean_obs = obs_discrep_mean;
        add_nug = @(X) X+eye(size(X))*1e-4;
        obs_discrep_vals = res.settings.obs_y - ...
            mean_sim(res.settings.obs_x(:,1),...
                repmat(mean(theta1_samps),...
                    size(res.settings.obs_x,1),1),...
                    res.settings.obs_x(:,2));
        prior_cov = @(rho,xo,t2o,xp,t2p,lambda) ... 
            gp_cov(rho,[xo ones(size(xo,1),1).*t2o],...
            [xp ones(size(xp,1),1).*t2p],lambda,false);
        updated_mean = @(y,xo,t2o,xp,t2p,rho,lambda,sigma2) ...
            mean_obs(xp,t2p) + ...
            (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
                add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
                    sigma2*eye(size(xo,1))))*...
            (y - mean_obs(xo,t2o)) ;
        updated_cov = @(xo,t2o,xp,t2p,rho,lambda,sigma2) ...
            add_nug(prior_cov(rho,xp,t2p,xp,t2p,lambda)) - ...
            (add_nug(prior_cov(rho,xp,t2p,xo,t2o,lambda)) / ...
                add_nug(prior_cov(rho,xo,t2o,xo,t2o,lambda) + ...
                    sigma2*eye(size(xo,1)))) * ...
            add_nug(prior_cov(rho,xo,t2o,xp,t2p,lambda)) ;
        additional_discrep_mean = @(t2) updated_mean(...
            obs_discrep_vals,...
            res.settings.obs_x(:,1),...
            res.settings.obs_x(:,2),...
            CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
            mean(res.obs_rho(burn_in:end,:)),...
            mean(res.obs_lambda(burn_in:end,:)),...
            mean(res.obs_var(burn_in:end)));
        additional_discrep_cov = @(x,t2) updated_cov(...
            res.settings.obs_x(:,1),...
            res.settings.obs_x(:,2),...
            CTO_obs_x_01, repmat(t2,size(CTO_obs_x_01,1),1),...
            mean(res.obs_rho(burn_in:end,:)),...
            mean(res.obs_lambda(burn_in:end,:)),...
            mean(res.obs_var(burn_in:end))); 

        % Get settings
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            CTO_obs_x,CTO_obs_t2,CTO_obs_y,CTO_des_x,CTO_des_y,...
            'min_x',xmin,'range_x',xrange,...
            'min_t1',t2min,'range_t1',t2range,'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',des_discrep,...
            'obs_discrep_mean',CTO_obs_discrep_mean,...
            'des_discrep',false,...
            'emulator_use',false,...
            'modular',false,...
            'EmulatorMean',CTO_mean_sim,...
            'obs_var',des_var,...
            'additional_discrep_mean',additional_discrep_mean,...
            'additional_discrep_cov',additional_discrep_cov,...
            'obs_rho_beta_params',des_rho_beta_params,...
            'obs_lambda_gam_params',des_lambda_gam_params,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',true,...
            'CTO',true);

        % Perform CTO calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        % Extract results
        ctod_results.theta2(:,ii) = res.theta1 ;
        ctod_results.theta2_hat(ii) = mean(res.theta1(burn_in:end));
        ctod_results.des_rho(:,ii) = res.obs_rho ;
        ctod_results.des_rho_hat(ii) = mean(res.obs_rho(burn_in:end));
        ctod_results.des_lambda(:,ii) = res.obs_lambda ;
        ctod_results.des_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));        
        ctod_results.theta2_var(ii) = var(res.theta1(burn_in:end)) ;
        ctod_results.des_rho_var(ii) = ...
            var(res.obs_rho(burn_in:end)) ;
        ctod_results.des_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        ctod_results.des_var(:,ii) = res.obs_var;
        ctod_results.des_var_hat(ii) = mean(res.obs_var(burn_in:end));
        ctod_results.settings = settings;

        
        all_results.dcto_results = dcto_results;
        all_results.khtc_results = khtc_results;
        all_results.ctod_results = ctod_results;
        
        save('temp','all_results');
        
        % Close windows every now and again
        if mod(ii,5) == 0, close all ; end
        
        % Update location in the loop
        fprintf('####################################\n');
        fprintf('COMPLETED STEP %g/%g OF DISCREP %g\n',ii,m,discrep);
        fprintf('####################################\n');
        
    end
    
    % Get optimal theta2 and corresponding theta1 value
    % Get and plot true theta2
    fmfn =@(z) dual_calib_example_fn(...
        .75,xmin,xrange,theta1,t1min,t1range,...
        z,t2min,t2range,0,1,discrep,false);
    % Try two different start locations
    [theta2_1, fval_1] = fmincon(fmfn,t2min+t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [theta2_2, fval_2] = fmincon(fmfn,t2min+3*t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [~,idx] = min([fval_1 fval_2]);
    theta2s = [theta2_1 theta2_2] ; theta2 = theta2s(idx) ; 
    
    % Add true parameter values, discrepancy version to results
    dcto_results.discrepancy = discrep;
    khtc_results.discrepancy = discrep;
    ctod_results.discrepancy = discrep;
    dcto_results.true_theta1 = theta1;
    khtc_results.true_theta1 = theta1;
    ctod_results.true_theta1 = theta1;
    dcto_results.true_theta2 = theta2;
    khtc_results.true_theta2 = theta2;
    ctod_results.true_theta2 = theta2;
    
    % Add results to full set of results
    results{discrep+1,1} = dcto_results;
    results{discrep+1,2} = khtc_results;
    results{discrep+1,3} = ctod_results;
    
    % Save
    locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-11-20_DCTO_vs_KOHCTO_results']);
    save(locstr,'results');
    
    % Close all those windows
    close all ;
    
    % Update location in the loop
    fprintf('####################################\n');
    fprintf('COMPLETED DISCREPANCY VERSION %g\n',discrep);
    fprintf('####################################\n');
    
end

%% Gather results for SDoE DCTO w/ estimation of obs_var, des_var
clc ; clearvars -except dpath ; close all ; 

for discrep = 0:6 % For each version of the discrepancy

    % Turn on observation discrepancy in the model iff discrep ~=0
    if discrep == 0
        obs_discrep = false;
    else 
        obs_discrep = true;
    end
    
    % Number of chains
    m = 2;

    for ii = 1:m
        % Set initial and final observation set sizes
        obs_initial_size = 5;
        obs_final_size = 20;

        % Settings
        modular = true;
        informative_targets = false;
        des_discrep = false;
        obs_discrep_use_MLEs = false;
        obs_var = 0.05; % observation error
        des_x_size = 15;
        doplot = true;
        verbose = false;
        obs_var_est = false;

        % Set number of draws, burn_in for each mcmc:
        M = 8e3; b = .5 ; burn_in = M*b;

        % Define inputs mins and ranges 
        xmin = .5;
        xrange = .5;
        t1min = 1.5;
        t1range = 3;
        t2min = 0;
        t2range = 5;

        % Observation and target discrepancy covariance hyperparameters
        obs_rho_beta_params = [8,1];
        obs_lambda_gam_params = [8,4];
        des_rho_beta_params = [2,.4];
        des_lambda_gam_params = [40,40];

        % True theta1 function
        high_theta1 = 2.25 ; low_theta1 = 1.5 ; 
%         t1_fn = @(x) high_theta1 - ...
%             (high_theta1-low_theta1) * ...
%             exp(40*((x-t2min)/t2range)-20)./...
%             (1+exp(40*((x-t2min)/t2range)-20));
        t1_fn = @(x) high_theta1 - ...
            (high_theta1-low_theta1) * ...
            exp(40*((x-t2min)/t2range)-20)./...
            (1+exp(40*((x-t2min)/t2range)-20));


        % Define comp. model & truth (version w/ nonstandardized output)
        model_fn_ns = @(x,t1,t2) dual_calib_example_fn(...
            x,xmin,xrange,t1,t1min,t1range,t2,t2min,t2range,...
            0,1,0,true) ; 
        true_phenomenon_ns = @(x,t2) dual_calib_example_fn(...
            x,xmin,xrange,...
            (t1_fn(t2*t2range+t2min)-t1min)./t1range,t1min,t1range,...
            t2,t2min,t2range,0,1,discrep,true);
            

        % Get mean and std using comp. model, define standardized version
        X=lhsdesign(1000,3); 
        Y=model_fn_ns(X(:,1),X(:,2),X(:,3));
        mean_y = mean(Y) ; std_y = std(Y) ;
        model_fn = @(x,t1,t2) (model_fn_ns(x,t1,t2) - mean_y)./std_y;
        true_phenomenon = ...
            @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

        % Get initial design
        X = lhsdesign(obs_initial_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless+ ...
            randn(obs_initial_size,1) * sqrt(obs_var) * std_y;

        % Now set desired observations
        des_x = linspace(0,1,des_x_size)' * xrange + xmin;
        if informative_targets    
            % Find true minimum at des_x
            fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
                t1_fn(t2),t1min,t1range,...
                t2,t2min,t2range,0,1,discrep,false);
            theta2 = ...
                fmincon(fn,rand*t2range + t2min,...
                [],[],[],[],t2min,t2min+t2range);
                true_min = dual_calib_example_fn(des_x,xmin,xrange,...
                    t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
                    theta2,t2min,t2range,0,1,discrep,false);
                des_y = true_min - 2*sqrt(obs_var) ; 
        else
            des_y = 0 * ones(size(des_x,1),1);
            % Get est of target error variance
            des_y_std = (des_y - mean_y)/std_y;
            des_var = (min(des_y_std)-min((Y-mean_y)/std_y))^2;
        end

        % Since we are not using emulator, empty out simulator observations
        sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

        % And get the average discrepancy value to use as the mean
        int_fn =@(x,t) true_phenomenon(x,t) - ...
            model_fn(x,(t1_fn(t*t2range+t2min)-t1min)./t1range,t) ;
        avg_disc = integral2(int_fn,0,1,0,1) ; 
        fprintf('Average observation discrepancy: %g\n',avg_disc);
        obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

        % Emulator mean
        mean_sim = @(a,b,c) model_fn(a,b,c) ; 

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,...
            des_x,des_y,'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_final_size',obs_final_size,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'true_phenomenon',true_phenomenon,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est,...
            'des_var_est',true,...
            'true_obs_var',obs_var);

        % Perform dual calibration
        % We need a loop because sometimes an ill-conditioned matrix early
        % in the burn-in makes the whole calibration fail.
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        
        sdoe_results.theta1(:,ii) = res.theta1 ;
        sdoe_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        sdoe_results.theta2(:,ii) = res.theta2 ;
        sdoe_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        sdoe_results.obs_rho(:,:,ii) = res.obs_rho ;
        sdoe_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        sdoe_results.obs_lambda(:,ii) = res.obs_lambda ;
        sdoe_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        sdoe_results.des_rho(:,ii) = res.des_rho ;
        sdoe_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        sdoe_results.des_lambda(:,ii) = res.des_lambda ;
        sdoe_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        sdoe_results.obs_var(:,ii) = res.obs_var ;
        sdoe_results.obs_var_hat(ii) = mean(res.obs_var(burn_in:end));
        sdoe_results.des_var(:,ii) = res.des_var ;
        sdoe_results.des_var_hat(ii) = mean(res.des_var(burn_in:end));
        sdoe_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        sdoe_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        sdoe_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        sdoe_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        sdoe_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        sdoe_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        sdoe_results.settings{ii} = res.settings ; 
        
        %%% Now do DCTO with all observations made ahead of time
        % Get observations
        X = lhsdesign(obs_final_size,2);
        obs_x = X(:,1) * xrange + xmin ; 
        obs_t2 = X(:,2) * t2range + t2min ;

        % Make a col vector based on true theta1
        obs_t1 = t1_fn(obs_t2);

        % Get "real" observations without noise but with discrepancy
        obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
            (obs_t2-t2min)./t2range);

        % Now add N(0,obs_var) noise to STANDARDIZED observations
        obs_y = obs_y_noiseless + ...
            randn(obs_final_size,1) * sqrt(obs_var) * std_y;

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
            obs_x,obs_t2,obs_y,des_x,des_y,...
            'min_x',xmin,'range_x',xrange,...
            'min_t1',t1min,'range_t1',t1range,...
            'min_t2',t2min,'range_t2',t2range,...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',obs_discrep,...
            'obs_discrep_mean',obs_discrep_mean,...
            'emulator_use',false,...
            'EmulatorMean',mean_sim,...
            'modular',modular,...
            'des_discrep',des_discrep,...
            'des_var',des_var,...
            'obs_var',obs_var,...
            'obs_rho_beta_params',obs_rho_beta_params,...
            'obs_lambda_gam_params',obs_lambda_gam_params,...
            'des_rho_beta_params',des_rho_beta_params,...
            'des_lambda_gam_params',des_lambda_gam_params,...
            'obs_discrep_use_MLEs',obs_discrep_use_MLEs,...
            'obs_final_size',obs_final_size,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est,...
            'des_var_est',true);

        % Perform dual calibration
        count = 0 ; err_count = 0 ; 
        while count == err_count
            try
                res = MCMC_dual_calib(settings);
            catch ME
                warning('Warning: calibration failed. Retrying...');
                err_count = err_count + 1;
            end
            count = count + 1;
            if count >= 10, rethrow(ME) ; end
        end
        
        pdoe_results.theta1(:,ii) = res.theta1 ;
        pdoe_results.theta1_hat(ii) = mean(res.theta1(burn_in:end)) ;
        pdoe_results.theta2(:,ii) = res.theta2 ;
        pdoe_results.theta2_hat(ii) = mean(res.theta2(burn_in:end));
        pdoe_results.obs_rho(:,:,ii) = res.obs_rho ;
        pdoe_results.obs_rho_hat(ii,:) = ...
            mean(res.obs_rho(burn_in:end,:));
        pdoe_results.obs_lambda(:,ii) = res.obs_lambda ;
        pdoe_results.obs_lambda_hat(ii) = ...
            mean(res.obs_lambda(burn_in:end));
        pdoe_results.des_rho(:,ii) = res.des_rho ;
        pdoe_results.des_rho_hat(ii) =mean(res.des_rho(burn_in:end));
        pdoe_results.des_lambda(:,ii) = res.des_lambda ;
        pdoe_results.des_lambda_hat(ii) = ...
            mean(res.des_lambda(burn_in:end));
        pdoe_results.obs_var(:,ii) = res.obs_var ;
        pdoe_results.obs_var_hat(ii) = mean(res.obs_var(burn_in:end));
        pdoe_results.des_var(:,ii) = res.des_var ;
        pdoe_results.des_var_hat(ii) = mean(res.des_var(burn_in:end));
        pdoe_results.theta1_var(ii) = var(res.theta1(burn_in:end)) ;
        pdoe_results.theta2_var(ii) = var(res.theta2(burn_in:end)) ;
        pdoe_results.obs_rho_var(ii,:) = ...
            var(res.obs_rho(burn_in:end,:)) ;
        pdoe_results.obs_lambda_var(ii) = ...
            var(res.obs_lambda(burn_in:end)) ;
        pdoe_results.des_rho_var(ii) = ...
            var(res.des_rho(burn_in:end)) ;
        pdoe_results.des_lambda_var(ii) = ...
            var(res.des_lambda(burn_in:end)) ;
        pdoe_results.settings{ii} = res.settings ; 
        
        both_results.sdoe_results = sdoe_results;
        both_results.pdoe_results = pdoe_results;
        
        save('temp','both_results');
        
        % Close windows every now and again
        if mod(ii,10) == 0, close all ; end
        
        % Update location in the loop
        fprintf('\n####################################\n');
        fprintf('COMPLETED STEP %g/%g OF DISCREP %g\n',ii,m,discrep);
        fprintf('####################################\n');
        
        
    end
    
    % Get optimal theta2 and corresponding theta1 value
    % Get and plot true theta2
    fmfn =@(z) dual_calib_example_fn(...
        .75,xmin,xrange,t1_fn(z),t1min,t1range,...
        z,t2min,t2range,0,1,discrep,false);
    % Try two different start locations
    [theta2_1, fval_1] = fmincon(fmfn,t2min+t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [theta2_2, fval_2] = fmincon(fmfn,t2min+3*t2range/4,...
        [],[],[],[],t2min,t2min+.5*t2range);
    [~,idx] = min([fval_1 fval_2]);
    theta2s = [theta2_1 theta2_2] ; theta2 = theta2s(idx) ; 
    theta1 = t1_fn(theta2);
    
    % Add true parameter values, discrepancy version to results
    sdoe_results.discrepancy = discrep;
    pdoe_results.discrepancy = discrep;
    sdoe_results.true_theta1 = theta1;
    pdoe_results.true_theta1 = theta1;
    sdoe_results.true_theta2 = theta2;
    pdoe_results.true_theta2 = theta2;
    
    % Add results to full set of results
    results{discrep+1,1} = sdoe_results;
    results{discrep+1,2} = pdoe_results;
    
    % Save
    locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-10-31_SDOE_results_desvarest_nobs' ...
    int2str(obs_initial_size) '-'...
    int2str(obs_final_size)]);
%     save(locstr,'results');
    
    % Close all those windows
    close all ;
    
    % Update location in the loop
    fprintf('####################################\n');
    fprintf('\nCOMPLETED DISCREPANCY VERSION %g\n',discrep);
    fprintf('####################################\n');
    
end
