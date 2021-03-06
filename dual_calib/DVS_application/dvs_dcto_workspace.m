% Dynamic Vibration System DCTO workspace
% Using data provided by Leslie Xu for DCTO
% 2020-09-11
%% Set path string and add paths
clc; clear all; close all;

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'dual_calib']);
addpath([dpath,'dual_calib\DVS_application']);
addpath([dpath,'dual_calib\DVS_application\data']);
addpath([dpath,'Example']);

% Change dir
cd(dpath);

% Let user know
fprintf('\nDone.\n')

%% Load the data
clc ; clearvars -except dpath ; close all ;

% Load the data
filename = [dpath,'dual_calib\DVS_application\data\FEM_0921_LeslieXu.xlsx'];
raw_dat = readtable(filename);

sim_dat = raw_dat(:,2:5);
exp_dat = raw_dat(1:12,10:(end-1));

% Save the simulation data and experimental data separately
sim_filename = [dpath,'dual_calib\DVS_application\data\simulation_data'];
exp_filename = [dpath,'dual_calib\DVS_application\data\experiment_data'];
% save(sim_filename,'sim_dat');
% save(exp_filename,'exp_dat');


%% Split data into separate variables and gather minima and ranges
clc ; clearvars -except dpath ; close all ;

% Load the raw data
sim_filename = [dpath,'dual_calib\DVS_application\data\simulation_data'];
exp_filename = [dpath,'dual_calib\DVS_application\data\experiment_data'];
load(sim_filename);
load(exp_filename);

% Convert tables to matrices
sim_dat = sim_dat{:,:};
exp_dat = exp_dat{:,:};

% Split up the data
sim_x = sim_dat(:,1);
sim_t1 = sim_dat(:,2);
sim_t2 = sim_dat(:,3);
sim_y = sim_dat(:,4);
obs_x = exp_dat(:,1);
obs_t2 = exp_dat(:,2);
obs_y = exp_dat(:,3);

% Get minima and ranges of inputs
x_min = min([sim_x;obs_x]);
x_range = range([sim_x;obs_x]);
t1_min = min(sim_t1);
t1_range = range(sim_t1);
t2_min = min([sim_t2;obs_t2]);
t2_range = range([sim_t2;obs_t2]);

% Get mean and std of outputs
y_mean = mean([sim_y;obs_y]);
y_std = std([sim_y;obs_y]);

% % Now set desired observations
des_x = linspace(0,1,8)' * x_range + x_min;
des_y = zeros(size(des_x,1),1);

% Now package everything up and save it
clear exp_dat sim_dat exp_filename sim_filename dpath
% save(['C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\',...
%     'dual_calib\DVS_application\data',...
%     '\2020-09-22_dvs_data_and_scaling_params']);
dpath = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\';

%% Use polynomial regression to get mean function for GP prior
clc ; clearvars -except dpath ; close all ;

% Load the raw data
loadfile = [dpath, 'dual_calib\DVS_application\data',...
    '\2020-09-22_dvs_data_and_scaling_params'];
load(loadfile);

% Get the mean function, mean_sim
mod_terms = []; 
coefs=[];
sim_x_01 = (sim_x - x_min)./x_range;
sim_t1_01 = (sim_t1 - t1_min)./t1_range;
sim_t2_01 = (sim_t2 - t2_min)./t2_range;
sim_y_std = (sim_y - y_mean)./y_std;
for ii=1:size(sim_y,2)
    mod = polyfitn([sim_x_01 sim_t1_01 sim_t2_01],sim_y_std(:,ii),2);
    mod_terms(:,:,ii) = mod.ModelTerms;
    coefs(:,ii) = mod.Coefficients;
end
mean_sim = @(xt_01) reshape(cell2mat(...
        arrayfun(@(row_idx) ...
            sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
                size(mod_terms,1),size(mod_terms,3)),1),...
            (1:size(xt_01,1)),...
            'UniformOutput',false)'),...
    size(xt_01,1),size(sim_y,2)) ;

clearvars -except loadfile mean_sim dpath;
load(loadfile)
clear loadfile
% save([dpath, 'dual_calib\DVS_application\data',...
%     '\2020-09-22_dvs_data_and_scaling_params']);



%% Get MLEs for GP surrogate
clc ; clearvars -except dpath ; close all ; 

% Load the raw data
loadfile = [dpath, 'dual_calib\DVS_application\data',...
    '\2020-09-22_dvs_data_and_scaling_params'];
load(loadfile);

% Normalize the simulation inputs, standardize the outputs
x = (sim_x - x_min) ./ x_range ;
t1 = (sim_t1 - t1_min) ./ t1_range;
t2 = (sim_t2 - t2_min) ./ t2_range;
y = (sim_y - y_mean) ./ y_std;

% Define function for minimization
f = @(rl) ...
    -logmvnpdf(y',...
    mean_sim([x t1 t2])',...
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

% Result w 0 mean: 
%   rho vals 0.609291544819824, 0.000019344581551, 0.125055951050202,
%   lambda val 0.710599411968917
% With deg-2 polynomial regression mean:
%   rho vals 0.358985186722316, 0.000000000038066, 0.026949566520274,
%   lambda val 5.843058403363464
inp


%% Perform dual calibration with emulator and discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
loadfile = [dpath, 'dual_calib\DVS_application\data',...
    '\2020-09-22_dvs_data_and_scaling_params'];
load(loadfile);

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,...
    des_x,des_y,...
    'min_x',x_min,'range_x',x_range,...
    'min_t1',t1_min,'range_t1',t1_range,...
    'min_t2',t2_min,'range_t2',t2_range,...
    'mean_y',y_mean,'std_y',y_std,...
    'M',1e4,...
    'burn_in',.5,...
    'EmulatorMean',mean_sim,...
    'obs_discrep',true,...
    'modular',true,...
    'obs_rho_beta_params',[10 .1],...
    'obs_lambda_gam_params',[5,5]);

% Perform calibration
tic;
results = MCMC_dual_calib(settings);
time_elapsed = toc;
results.time_elapsed = time_elapsed;

% Save
saveloc = [dpath, 'dual_calib\DVS_application\data',...
    '\2020-09-30_dvs_dcto_results'];
% save(saveloc,'results');


%% Gather NSGA-II results for DVS
clc ; clearvars -except dpath ; close all ;

% Load the cost grid results
% locstr= ['C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\'...
%     'stored_data\2020-04-26_CTO_costgrid_size30'];
locstr = [dpath, 'dual_calib\DVS_application\data',...
    '\2020-09-30_dvs_dcto_results'];
load(locstr);

% Prepare elements needed for bounds and objective function
sim_t2 = results.settings.sim_t2 .* results.settings.range_t2 +...
    results.settings.min_t2;

% To use a grid of three x values, uncomment this:
x = linspace(0,1,3)' .* results.settings.range_x + results.settings.min_x;

% To use same des_x as DCTO, uncomment this:
% x = results.settings.des_x .* results.settings.range_x + ...
%     results.settings.min_x;

% Store size of x
sizex = size(x,1);

% Get posterior point estimate of calibration input
t1 = mean(results.theta1);

% Set NSGA-II options
options = nsgaopt();                    % create default options structure
options.popsize = 50;                   % population size
options.maxGen  = 100;                  % max generation

options.numObj = sizex;                     % number of objectives
options.numVar = size(results.theta2,2); % number of design variables
options.numCons = 0;                    % number of constraints
options.lb = results.settings.min_t2;               % lower bound of x
options.ub = results.settings.min_t2 + results.settings.range_t2;%up bd x
options.objfun = @(t2) emulator_mean(results,x,ones(sizex,1).*[t1 t2]);% objective function handle

options.plotInterval = 10;              % large interval for efficiency
options.outputInterval = 10;

tic;
nsga_result = nsga2(options);
elapsed_time = toc;
nsga_result.elapsed_time = elapsed_time;


% For convenience, add top-level array to result with final generation
% design settings and objective values.
final_theta2 = nan(options.popsize,options.numVar);
final_obj = nan(options.popsize,options.numObj);

for ii=1:options.popsize
    final_theta2(ii,:) = nsga_result.pops(options.maxGen,ii).var;
    final_obj(ii,:) = nsga_result.pops(options.maxGen,ii).obj;
end

nsga_result.final_theta2 = final_theta2;
nsga_result.final_obj = final_obj;

results.nsga_result=nsga_result;
    
% Output the results:
unique(nsga_result.final_theta2)
unique(nsga_result.final_obj,'rows' )

% Save NSGA-II result
% save(locstr,'results');
