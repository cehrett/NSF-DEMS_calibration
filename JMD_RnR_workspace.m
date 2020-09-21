% Results gathered for use in the CTO paper sent as revised version to JMD.
% This file was not used for results presented in the draft that was
% submitted to Technometrics. I created this file in order to gather
% results when revising the paper after the initial reviewer comments. 
% I began working on revisions using the CTO_results_workspace.m file so it
% may be that not all new results are gathered here, some may be there.

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
addpath([dpath,'Example']);
addpath([dpath,'Example\Ex_results']);
addpath([dpath,'dual_calib']);
fprintf('Ready.\n');

% Change dir
cd(dpath);

%% ZDT Get mins, ranges, means, stds from objective fn
clc ; clearvars -except dpath ; close all ;

% Utopia point: [0 0]
% Dystopia point: [1+.5^Len(x), 1+.5^Len(x) ]

X=rand(1e6,5);
Y = TP_ZDT1_objfun(X);
mean_y = mean(Y);
std_y = std(Y);

% Get the mean and std of the distance from the mean
dists = sqrt(sum((Y-mean_y).^2,2));
mean_dist = mean(dists);
std_dist = std(dists);

% Get mean and s.d. of the distance from the mean on standardized scale
Y_std = (Y-mean_y)./std_y ;
dists_std = sqrt(sum(Y_std.^2,2));
mean_dist_std = mean(dists_std);
std_dist_std = std(dists_std);


%% ZDT Gather PCTO results, estimating var, no emulator
clc ; clearvars -except dpath ; close all ;

% Settings
obs_var = 5000; % initial guess of observation error var
num_target_obs = 1;
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 2e4; b = .5 ;

% Define inputs mins and ranges 
xmin = [];
xrange = [];
tmin = [0 0 0 0 0];
trange = [1 1 1 1 1];

% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(xt) ...
    TP_ZDT1_objfun(xt .* [xrange trange] + [xmin tmin]);

% Objective function mean and std were found via brute force (elsewhere)
mean_y = [0.499976319936033,3.947381285207361] ; 
std_y = [0.288667391125751,1.240743941520046] ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% Get target outcomes
obs_x = zeros(num_target_obs,size(xmin,2));
% The following line sets the target to be the utopia point of the system
obs_y = repmat([0 0],num_target_obs,1); 
% Other possibilities not used here:
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% obs_y = repmat([0.7130 0.7144 17.9220],size(obs_x,1),1);

% Emulator mean: set to true function (standardized scale)
mean_sim = @(x,t1,t2) model_fn([x t1 t2]);
% If we were using an emulator, we would set:
% mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings([],[],[],[],...
    obs_x,[],obs_y,...
    [],zeros(0,3),'min_x',obs_x,'range_x',obs_x,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',false,...
    'CTO',true);

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

% save([dpath,'Example\Ex_results\'...
%     '2020-04-11_PCTO_TP_ZDT1'],...
%     'res');


%% ZDT Gather results without PCTO, estimating var, no emulator
clc ; clearvars -except dpath ; close all ;

% Settings
obs_var = 0.5; % initial guess of observation error var
num_target_obs = 1;
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 2e4; b = .5 ;

% Define inputs mins and ranges 
xmin = [];
xrange = [];
tmin = [0 0 0 0 0];
trange = [1 1 1 1 1];

% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(xt) ...
    TP_ZDT1_objfun(xt .* [xrange trange] + [xmin tmin]);

% Objective function mean and std were found via brute force (elsewhere)
mean_y = [0.499976319936033,3.947381285207361] ; 
std_y = [0.288667391125751,1.240743941520046] ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% Get target outcomes
obs_x = zeros(num_target_obs,size(xmin,2));
% The following line sets the target using target selection code elsewhere
% obs_y = repmat([0.25   -6],num_target_obs,1); % 5.222 SD from feasible
obs_y = repmat([0.3984   -2.1501],num_target_obs,1); % 2 SD from feasible
% Other possibilities not used here:
% obs_y = repmat([0 0],num_target_obs,1); % Utopia pt
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% obs_y = repmat([0.7130 0.7144 17.9220],size(obs_x,1),1);

% Emulator mean: set to true function (standardized scale)
mean_sim = @(xt) model_fn(xt);
% If we were using an emulator, we would set:
% mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings([],[],[],[],...
    obs_x,[],obs_y,...
    [],zeros(0,2),'min_x',obs_x,'range_x',obs_x,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'obs_var_same',true,...
    'CTO',true);

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

% locstr = [dpath,'Example\Ex_results\'...
%     '2020_04_19_CTO_noemulator_TP_ZDT1_5p222SD'];
locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_TP_ZDT1_2SD'];
% save(locstr,'res');

%% ZDT Identify nearby target outcome (without PCTO)
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = [];
xrange = [];
tmin = [0 0 0 0 0];
trange = [1 1 1 1 1];

% Objective function mean and std were found via brute force (elsewhere)
mean_y = [0.499976319936033,3.947381285207361] ; 
std_y = [0.288667391125751,1.240743941520046] ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% The following line sets the target to be the utopia point of the system,
% which was located using fmincon on each objective function
target_original = [0.25   -6 ]; 
target_original_std = (target_original - mean_y)./std_y;

% We will find the output that minimizes the distance to the utopia point,
% using fmincon. All outputs on standardized scale.
minfn = @(x) sum(((TP_ZDT1_objfun(x)-mean_y)./std_y-target_original_std).^2);
[X,~] = fmincon(minfn,[.5 .5 .5 .5 .5],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);

% Now we can get the feasible target
target_feasible = TP_ZDT1_objfun(X);
target_feasible_std = (target_feasible - mean_y)./std_y;

% Move to a target 1 s.d. away from the feasible target, in the direction
% of the original target
direction_vector = target_feasible_std - target_original_std;
direction_vector_normd = direction_vector / norm(direction_vector);

% This is the distance we want the new obs to be from the feasible region:
std_dists =  2;

% New target (on standarized scale and original scale)
target_new_std = target_feasible_std - std_dists * direction_vector_normd ;
target_new = target_new_std .* std_y + mean_y

% Also get the Euclidean distances of the standardized targets from the
% feasible point
target_original_distance = ...
    sqrt(sum((target_original_std - target_feasible_std).^2));
target_new_distance = ...
    sqrt(sum((target_new_std - target_feasible_std).^2));



%% ZDT Get optimal results using fmincon
clc ; clearvars -except dpath ; close all

tar_obj = [ .25 -6 ];

mean_y = [0.499976319936033,3.947381285207361] ; 
std_y = [0.288667391125751,1.240743941520046] ;

fun = @(x) sqrt(sum(((TP_ZDT1_objfun(x)-mean_y)./std_y - (tar_obj-mean_y)./std_y).^2,2));
x0 = rand(1,5);
opt_des = fmincon(fun,x0,[],[],[],[],[0 0 0 0 0],[1 1 1 1 1]);

opt_des
opt_obj = fun(opt_des)

x1 = rand(1e6,5);
dists = fun(x1);
sum(dists < opt_obj)
fun([.5 0.1 0.1 0.1 0.1 ])


%% WTA Gather NSGA-II results for the wind turbine application
clc ; clearvars -except dpath ; close all ;

% Load the cost grid results
locstr = [dpath,'stored_data\'...
    '2019-11-05_CTO_costgrid'];
load(locstr);

options = nsgaopt();                    % create default options structure
options.popsize = 50;                   % populaion size
options.maxGen  = 500;                  % max generation

options.numObj = 2;                     % number of objectives
options.numVar = 2;                    % number of design variables
options.numCons = 0;                    % number of constraints
options.lb = [.1 10];               % lower bound of x
options.ub = [.6 25];                % upper bound of x
options.objfun = @(x) emulator_mean(results{1},x);% objective function handle

options.plotInterval = 10;              % large interval for efficiency
options.outputInterval = 10;

tic;
result = nsga2(options);
elapsed_time = toc;
result.elapsed_time = elapsed_time;

% For convenience, add top-level array to result with final generation
% design settings and objective values.
final_theta = nan(options.popsize,options.numVar);
final_obj = nan(options.popsize,options.numObj);

for ii=1:options.popsize
    final_theta(ii,:) = result.pops(options.maxGen,ii).var;
    final_obj(ii,:) = result.pops(options.maxGen,ii).obj;
end

result.final_theta = final_theta;
result.final_obj = final_obj;

% Save NSGA-II result
locstr = [dpath,'stored_data\'...
    '2020-04-23_NSGA2_results'];
% save(locstr,'result');

%% WTA_alt Perform CTO on bivariate objective over cost grid (for PF est)
clc ; clearvars -except dpath ; close all ;

% Store all results here:
all_results = {};

%%% Load raw data
locstr = [dpath,'stored_data\'...
    '2020-04-23_new_WTA_cases_size100'];
load(locstr);

for jj=1:length(dats)

    sim_x = dats{jj}.sim_x;
    sim_t = dats{jj}.sim_t;
    sim_y = dats{jj}.sim_y; 


    % Define inputs mins and ranges 
    xmin = min(sim_x);
    xrange = range(sim_x);
    tmin = min(sim_t);
    trange = range(sim_t);

    % Settings
    obs_x_size  = 1; % Number of target observations
    doplot = true;
    verbose = true;

    % Set number of draws, burn_in for each mcmc:
    M = 1e4; b = .5 ;

    % Get output means and stds
    mean_y = mean(sim_y) ; std_y = std(sim_y) ;

    % Get target design
    obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;

    % Emulator mean
%     mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 
    % Evaluate with deg 2 poly regression
    mod_terms = [];
    coefs=[];
    sim_x_01 = (sim_x - xmin)./xrange;
    sim_t_01 = (sim_t - tmin)./trange;
    sim_y_std = (sim_y - mean_y)./std_y;
    for ii=1:size(sim_y,2)
        mod = polyfitn([sim_x_01 sim_t_01],sim_y_std(:,ii),2);
        mod_terms(:,:,ii) = mod.ModelTerms;
        coefs(:,ii) = mod.Coefficients;
    end
    mean_sim = @(xt_01) reshape(cell2mat(...
        arrayfun(@(row_idx) ...
        sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
        size(mod_terms,1),size(mod_terms,3)),1),(1:size(xt_01,1)),...
        'UniformOutput',false)'),size(xt_01,1),size(sim_y,2)) ;

    % Define grid of targets
    ngrid = 20 ; % Number of points in grid
    % Minimum and maximum possible cost:
    mincost = min(sim_y(:,2)); maxcost = max(sim_y(:,2)); 
    mindefl = min(sim_y(:,1)); maxdefl = max(sim_y(:,1)); 
    des_obs_grid = [zeros(ngrid/2,1) linspace(mincost,maxcost,ngrid/2)'; ...
        linspace(mindefl,maxdefl,ngrid/2)' zeros(ngrid/2,1)] ;
    obs_var_est_grid = [repmat([true false],ngrid/2,1);...
        repmat([false true],ngrid/2,1)];
    obs_var_grid = [repmat([1 0.00001],ngrid/2,1);...
        repmat([0.00001 1],ngrid/2,1)]; % initial guess of obs error var

    % Perform CTO for each point in the grid
    for ii = 1 : ngrid

        % Select the target design to be the point chosen in CTO (elsewhere)
        obs_y = repmat(des_obs_grid(ii,:),obs_x_size,1);

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
            obs_x,[],obs_y,...
            [],zeros(0,size(obs_y,2)),'min_x',xmin,'range_x',xrange,...
            'min_t1',tmin,'range_t1',trange,...
            'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',false,...
            'emulator_use',true,...
            'EmulatorMean',mean_sim,... 
            'EmulatorCovHypers','Default',...
            'des_discrep',false,...
            'obs_var',obs_var_grid(ii,:),...
            'obs_discrep_use_MLEs',false,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est_grid(ii,:),...
            'CTO',true);

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

        results{ii} = res;
        save temp

    end

    all_results{jj}=results;
    clear results
    
end

% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_CTO_double_grid_size',int2str(ngrid)];
% save(locstr,'all_results');

%% WTA_alt Get posterior predictive results for cost grid estimate of PF
clc ; clearvars -except dpath ; close all ;

% Load the results
% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
% locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-26_CTO_double_grid_size20';
% locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-29_3_CTO_grids_size20';
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-30_3_CTO_double_grid_size20';
% locstr = [dpath,'stored_data\'...
%     '2020-04-23_CTO_costgrid_new_cases_size100'];
load(locstr);

for hh=1:length(all_results)

    results = all_results{hh};

    % Define inputs mins and ranges 
    xmin = results{1}.settings.min_x;
    xrange = results{1}.settings.range_x;
    tmin = results{1}.settings.min_t1;
    trange = results{1}.settings.range_t1;

    % Build whatever matrices can be done outside of the main loop
    sim_xt = [results{1}.settings.sim_x results{1}.settings.sim_t1 ...
        results{1}.settings.sim_t2];
    eta = results{1}.settings.sim_y;
    mean_sim = results{1}.settings.mean_sim;
    prior_mean_simxt = mean_sim(sim_xt) ;
    cov_xt_xt = nan(size(eta,1),size(eta,1),2);
    burn_in = results{1}.settings.burn_in;
    for ii = 1 : 2
        emulator_rho = results{1}.settings.emulator_rho(:,ii);
        emulator_lambda = results{1}.settings.emulator_lambda(ii);
        cov_xt_xt(:,:,ii) = ...
            gp_cov(emulator_rho, sim_xt, sim_xt, ...
            emulator_lambda,false);
        cov_xt_xt(:,:,ii) = cov_xt_xt(:,:,ii) + ...
            1e-5*eye(size(cov_xt_xt(:,:,ii)));
    end

    for kk = 1 : length(results)

        % Get posterior samples
        theta = results{kk}.theta1;

        % Get set up for loop that will find posterior predictive dist
        pred_xt_all = ([.5*ones(size(theta,1),size(xmin,2)) theta] - ...
            [zeros(1,size(xmin,2)) tmin])./[ones(1,size(xmin,2)) trange];
        prior_mean_xt_all = mean_sim(pred_xt_all) ; 

        % Get posterior predictive distribution for each output
        % Only do maxsimul samples at once, to ease computational expense
        maxsimul = 1000;
        for jj = 1 : ceil(size(theta,1)/maxsimul)
            startpt = (jj-1) * maxsimul + 1 ;
            endpt = min(jj * maxsimul, size(theta,1));
            pred_xt = pred_xt_all(startpt:endpt,:);
            for ii=1:2
                emulator_rho = results{kk}.settings.emulator_rho(:,ii);
                emulator_lambda = results{kk}.settings.emulator_lambda(ii);
                cov_pxt_xt = ...
                    gp_cov(emulator_rho, pred_xt, sim_xt, ...
                    emulator_lambda,false);
                cov_pxt_pxt = ...
                    gp_cov(emulator_rho, pred_xt, pred_xt,...
                    emulator_lambda,false);
                pred_mean_xt_current(:,ii) = ...
                    prior_mean_xt_all(startpt:endpt,ii) + ...
                    cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
                    (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
                pred_cov_xt_current(:,:,ii) = ...
                    cov_pxt_pxt - cov_pxt_xt * ...
                        (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
                fprintf('\n%g of %g, %g of 3 Done\n',jj,...
                    ceil(size(theta,1)/maxsimul),ii);
                pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
                pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
                    pred_cov_xt_current(:,:,ii);
            end
        end
        % Put back on original scale
        mean_y = results{kk}.settings.mean_y; 
        std_y = results{kk}.settings.std_y;
        pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
        mean(pred_mean_xt_os)


        % Now resample from the relevant gaussian dists, for full uncertainty
        for ii = 1:2
            pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
            pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
                pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
        end
        % Put back on original scale
        pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
        mean(pred_samps_xt_os)


        %%%%%%%%%%% Make a figure
        % xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
        % ylims = [0   95 ; 0   120 ; 0    3.5000];
        figure('Position',[10 10 600 300],'Color','white');
        for ii = 1:2
            subplot(1,3,ii);
    %         histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
    %             'EdgeColor','none');
    %         hold on;
            histogram(pred_mean_xt_os(burn_in:end,ii),'Normalization','pdf',...
                'EdgeColor','none');
            xlabel(sprintf('y_%g',ii));
            set(gca,'Yticklabel',[]);
            set(gca,'Ytick',[]);
            drawnow;
        %     xlim(xlims(ii,:));
        %     ylim(ylims(ii,:));
        end

        % Save the posterior predictive output
        model_output.means = pred_mean_xt;
        model_output.vars = cell2mat(...
            arrayfun(@(b)diag(pred_cov_xt(:,:,b)),1:2,'UniformOutput',false));
        all_results{hh}{kk}.model_output = model_output;
%         save temp
    end
end

% save(locstr,'all_results');

%% WTA_alt Gather NSGA-II results for the new wind turbine applications
clc ; clearvars -except dpath ; close all ;

% Load the cost grid results
% locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-26_CTO_double_grid_size20';
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-29_3_CTO_grids_size20';
% locstr = [dpath,'stored_data\'...
%     '2020-04-23_CTO_costgrid_new_cases_size100'];
load(locstr);

for hh=1:length(all_results)
    
    results = all_results{hh};
    sim_t = results{1}.settings.sim_t1 .* results{1}.settings.range_t1 +...
        results{1}.settings.min_t1;

    options = nsgaopt();                    % create default options structure
    options.popsize = 50;                   % populaion size
    options.maxGen  = 500;                  % max generation

    options.numObj = 2;                     % number of objectives
    options.numVar = size(results{1}.theta1,2); % number of design variables
    options.numCons = 0;                    % number of constraints
    options.lb = min(sim_t);               % lower bound of x
    options.ub = max(sim_t);                % upper bound of x
    options.objfun = @(x) emulator_mean(results{1},x);% objective function handle

    options.plotInterval = 10;              % large interval for efficiency
    options.outputInterval = 10;

    tic;
    nsga_result = nsga2(options);
    elapsed_time = toc;
    nsga_result.elapsed_time = elapsed_time;

    % For convenience, add top-level array to result with final generation
    % design settings and objective values.
    final_theta = nan(options.popsize,options.numVar);
    final_obj = nan(options.popsize,options.numObj);

    for ii=1:options.popsize
        final_theta(ii,:) = nsga_result.pops(options.maxGen,ii).var;
        final_obj(ii,:) = nsga_result.pops(options.maxGen,ii).obj;
    end

    nsga_result.final_theta = final_theta;
    nsga_result.final_obj = final_obj;
    
    all_results{hh}{end+1}=nsga_result;
    
end

% Save NSGA-II result
% save(locstr,'all_results');

%% WTA_alts Load the new data
clc ; clearvars -except dpath ; close all 

% set random seed
rng(355);

% Load the data
dats = {};
dats{1}= xlsread('stored_data\2020-04-22_new_WTA_data',1);
dats{2}= xlsread('stored_data\2020-04-22_new_WTA_data',2);
dats{3}= xlsread('stored_data\2020-04-22_new_WTA_data2',1);
dats{4}= xlsread('stored_data\2020-04-22_new_WTA_data2',2);
dats{5}= xlsread('stored_data\2020-04-22_new_WTA_data2',3);
for ii=1:4
    raw_dat = dats{ii};
    k=100; raw_dat = raw_dat(sort(randsample(size(raw_dat,1),k)),:);
    sim_t = raw_dat(:,1:2);
    sim_y = raw_dat(:,[8,3]); 
    sim_x = zeros(size(sim_t,1),0);
    dats{ii}=struct('sim_t',sim_t,'sim_y',sim_y,'sim_x',sim_x);
    clear raw_dat
end
% Fifth case had different structure and must be handled separately
raw_dat = dats{5};
raw_dat = raw_dat(sort(randsample(size(raw_dat,1),k)),:);
sim_t = raw_dat(:,1:3);
sim_y = raw_dat(:,[9,4]); 
sim_x = zeros(size(sim_t,1),0);
dats{5}=struct('sim_t',sim_t,'sim_y',sim_y,'sim_x',sim_x);
clear raw_dat

% Save the processed data
locstr = [dpath,'stored_data\'...
    '2020-04-23_new_WTA_cases_size100'];
% save(locstr,'dats');

%% WTA Perform CTO using target chosen in PCTO
clc ; clearvars -except dpath ; close all ; 

%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
% Uncomment in order to reduce the size of the training set:
k=30; rng(355); raw_dat = raw_dat(sort(randsample(size(raw_dat,1),k)),:);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6); 
clear raw_dat

% Settings
obs_var = 0.5; % initial guess of observation error var
obs_x_size  = 3; % Number of target observations
doplot = true;
verbose = true;

% % Set emulator hyperparameter MLEs
% % These were found via fmincon.
% EmulatorCovHypers = [...
%     0.723919426870610   0.710378993535344   0.999999999999718;
%     0.978844923780470   0.972316694561107   0.998830926614838;
%     0.990609430210837   0.988186044400943   0.998556385291986;
%     0.017707123599447   0.026113351305598   0.000949737279367];
% Actually just estimate them at time of getting settings.
EmulatorCovHypers = 'Default';

% Set number of draws, burn_in for each mcmc:
M = 1e4; b = .5 ;

% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% Select the target design to be the point chosen in PCTO (elsewhere)
obs_y = repmat(...
    [0.742497187630386   0.089184738492403  71.112749103024015],...
    obs_x_size,1);

% Emulator mean
%     mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 
% Evaluate with deg 2 poly regression
mod_terms = [];
coefs=[];
sim_x_01 = (sim_x - xmin)./xrange;
sim_t_01 = (sim_t - tmin)./trange;
sim_y_std = (sim_y - mean_y)./std_y;
for ii=1:size(sim_y,2)
    mod = polyfitn([sim_x_01 sim_t_01],sim_y_std(:,ii),2);
    mod_terms(:,:,ii) = mod.ModelTerms;
    coefs(:,ii) = mod.Coefficients;
end
mean_sim = @(xt_01) reshape(cell2mat(...
    arrayfun(@(row_idx) ...
    sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
    size(mod_terms,1),size(mod_terms,3)),1),(1:size(xt_01,1)),...
    'UniformOutput',false)'),size(xt_01,1),size(sim_y,2)) ;


% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
    obs_x,[],obs_y,...
    [],zeros(0,3),'min_x',xmin,'range_x',xrange,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',true,...
    'EmulatorMean',mean_sim,...
    'EmulatorCovHypers',EmulatorCovHypers,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_var_same',true,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'CTO',true);

% Perform dual calibration
% We need a loop because sometimes an ill-conditioned matrix early
% in the burn-in makes the whole calibration fail.
count = 0 ; err_count = 0 ; 
while count == err_count
    try
        tic;
        res = MCMC_dual_calib(settings);
        elapsed_time = toc;
        res.elapsed_time = elapsed_time;
    catch ME
        warning('Warning: calibration failed. Retrying...');
        err_count = err_count + 1;
    end
    count = count + 1;
    if count >= 10, rethrow(ME) ; end
end

locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_CTO_size',int2str(k)];
% save(locstr,'res');

%% WTA Get posterior predictive distribution from CTO after PCTO
clc ; clearvars -except dpath ; close all ;

% Load the results
% If results already in memory can run this: clearvars -except dpath res
% locstr = [dpath,'stored_data\'...
%     '2020-04-25_CTO_size500'];
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-30_CTO_size30';
load(locstr);

% Define inputs mins and ranges 
xmin = res.settings.min_x;
xrange = res.settings.range_x;
tmin = res.settings.min_t1;
trange = res.settings.range_t1;

% Get posterior samples
burn_in = res.settings.burn_in;
t1 = res.theta1(:,1); t2 = res.theta1(:,2);

% Get set up for loop that will find posterior predictive dist
pred_xt_all = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
sim_xt = [res.settings.sim_x res.settings.sim_t1 res.settings.sim_t2];
prior_mean_xt_all = res.settings.mean_sim(pred_xt_all) ; 
eta = res.settings.sim_y;
% num_mod_obs = 8
% mod_obs_idx = randsample(1:504,num_mod_obs);
% sim_xt = sim_xt(mod_obs_idx,:);
% eta = eta(mod_obs_idx,:);
prior_mean_simxt = res.settings.mean_sim(sim_xt) ; 
cov_xt_xt = nan(size(sim_xt,1),size(sim_xt,1),3);
for ii=1:3
    cov_xt_xt(:,:,ii) = ...
        gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt(:,:,ii) = cov_xt_xt(:,:,ii) + ...
        1e-5*eye(size(cov_xt_xt(:,:,ii)));
end

% Get posterior predictive distribution for each output
% Only do maxsimul samples at once, to ease computational expense
maxsimul = 1000;
for jj = 1 : ceil(size(t1,1)/maxsimul)
    startpt = (jj-1) * maxsimul + 1 ;
    endpt = min(jj * maxsimul, size(t1,1));
    pred_xt = pred_xt_all(startpt:endpt,:);
    for ii=1:3
        cov_pxt_xt = ...
            gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
            res.settings.emulator_lambda(ii),false);
        cov_pxt_pxt = ...
            gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
            res.settings.emulator_lambda(ii),false);
        pred_mean_xt_current(:,ii) = ...
            prior_mean_xt_all(startpt:endpt,ii) + ...
            cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
            (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
        pred_cov_xt_current(:,:,ii) = ...
            cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
        fprintf('\n%g of %g, %g of 3 Done\n',jj,...
            ceil(size(t1,1)/maxsimul),ii);
        pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
        pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
            pred_cov_xt_current(:,:,ii);
    end
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(pred_mean_xt_os)

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(pred_samps_xt_os)


%%%%%%% Get prior predictive output 

% Number of samples
m = 1e4 ; 

% Get theta values
thetavals = rand(m,2) .* trange + tmin ; 
t1 = thetavals(:,1) ; t2 = thetavals(:,2); 

% Get set up for loop that will find prior predictive dist
pred_xt = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 

% Get posterior predictive distribution for each output
pred_mean_xt_prior = nan(m,3);
pred_cov_xt_prior = nan(m,m,3);
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt_prior(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt_prior(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
    fprintf('\n%g/3 done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
prior_mean_xt_os = pred_mean_xt_prior .* std_y + mean_y;
mean(prior_mean_xt_os)


% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt_prior(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt_prior,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt_prior(:,ii);
end
% Put back on original scale
prior_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(prior_samps_xt_os)


%%%%%%%%%%% Make a figure
% xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
% ylims = [0   95 ; 0   120 ; 0    3.5000];
fig1 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_mean_xt_os(5000:end,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_samps_xt_os(5000:end,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_without_GP_uncertainty' -m2

% Save the posterior predictive output
model_output.means = pred_mean_xt;
model_output.vars = cell2mat(...
    arrayfun(@(b)diag(pred_cov_xt(:,:,b)),1:3,'UniformOutput',false));
res.model_output = model_output;
% save(locstr,'res');

prior_model_output.theta = thetavals;
prior_model_output.means = pred_mean_xt_prior;
prior_model_output.vars = cell2mat(...
    arrayfun(...
        @(b)diag(pred_cov_xt_prior(:,:,b)),1:3,'UniformOutput',false));
locstr2 = [dpath,'stored_data\'...
    '2019-11-06_prior_predictive_distributions'];
% save(locstr2,'prior_model_output')

%% WTA Perform CTO on bivariate objective over cost grid (for PF est)
clc ; clearvars -except dpath ; close all ;


%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
% Uncomment in order to reduce the size of the training set:
k=30; rng(355); raw_dat = raw_dat(sort(randsample(size(raw_dat,1),k)),:);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,[4,6]); 
clear raw_dat


% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Settings
obs_var = [1 0.00001]; % initial guess of observation error var
obs_x_size  = 1; % Number of target observations
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 1.5e4; b = 2/3 ;

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;

% Emulator mean
%     mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 
% Evaluate with deg 2 poly regression
mod_terms = [];
coefs=[];
sim_x_01 = (sim_x - xmin)./xrange;
sim_t_01 = (sim_t - tmin)./trange;
sim_y_std = (sim_y - mean_y)./std_y;
for ii=1:size(sim_y,2)
    mod = polyfitn([sim_x_01 sim_t_01],sim_y_std(:,ii),2);
    mod_terms(:,:,ii) = mod.ModelTerms;
    coefs(:,ii) = mod.Coefficients;
end
mean_sim = @(xt_01) reshape(cell2mat(...
    arrayfun(@(row_idx) ...
    sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
    size(mod_terms,1),size(mod_terms,3)),1),(1:size(xt_01,1)),...
    'UniformOutput',false)'),size(xt_01,1),size(sim_y,2)) ;

% Define grid of targets
ngrid = 20 ; % Number of points in grid
% Minimum and maximum possible cost:
mincost = min(sim_y(:,2)); maxcost = max(sim_y(:,2)); 
des_obs_grid = [zeros(ngrid,1) linspace(mincost,maxcost,ngrid)'] ;

% Perform CTO for each point in the grid
for ii = 1 : ngrid

    % Select the target design to be the point chosen in CTO (elsewhere)
    obs_y = repmat(des_obs_grid(ii,:),obs_x_size,1);

    % Get settings for DCTO
    settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
        obs_x,[],obs_y,...
        [],zeros(0,size(obs_y,2)),'min_x',xmin,'range_x',xrange,...
        'min_t1',tmin,'range_t1',trange,...
        'min_t2',[],'range_t2',[],...
        'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
        'obs_discrep',false,...
        'emulator_use',true,...
        'EmulatorMean',mean_sim,... 
        'EmulatorCovHypers','Default',...
        'des_discrep',false,...
        'obs_var',obs_var,...
        'obs_discrep_use_MLEs',false,...
        'doplot',doplot,...
        'verbose',verbose,...
        'obs_var_est',[true false ],...
        'CTO',true);

    % Perform dual calibration
    % We need a loop because sometimes an ill-conditioned matrix early
    % in the burn-in makes the whole calibration fail.
    count = 0 ; err_count = 0 ; 
    while count == err_count
        try
            tic;
            res = MCMC_dual_calib(settings);
            elapsed_time = toc;
            res.elapsed_time = elapsed_time;
        catch ME
            warning('Warning: calibration failed. Retrying...');
            err_count = err_count + 1;
        end
        count = count + 1;
        if count >= 10, rethrow(ME) ; end
    end

    results{ii} = res;
    save temp

end


% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_CTO_costgrid_size',int2str(k)];
% save(locstr,'results');

%% WTA Get posterior predictive results for cost grid estimate of PF
clc ; clearvars -except dpath ; close all ;

% Load the results
% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
% locstr = [dpath,'stored_data\'...
%     '2020-04-23_CTO_costgrid_new_cases_size100'];
% locstr= ['C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\'...
%     'stored_data\2020-04-26_CTO_costgrid_size30'];
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-28_CTO_costgrid_size30';
load(locstr);


% Define inputs mins and ranges 
xmin = results{1}.settings.min_x;
xrange = results{1}.settings.range_x;
tmin = results{1}.settings.min_t1;
trange = results{1}.settings.range_t1;

% Build whatever matrices can be done outside of the main loop
sim_xt = [results{1}.settings.sim_x results{1}.settings.sim_t1 ...
    results{1}.settings.sim_t2];
eta = results{1}.settings.sim_y;
mean_sim = results{1}.settings.mean_sim;
prior_mean_simxt = mean_sim(sim_xt) ;
cov_xt_xt = nan(size(eta,1),size(eta,1),2);
burn_in = results{1}.settings.burn_in;
for ii = 1 : 2
    emulator_rho = results{1}.settings.emulator_rho(:,ii);
    emulator_lambda = results{1}.settings.emulator_lambda(ii);
    cov_xt_xt(:,:,ii) = ...
        gp_cov(emulator_rho, sim_xt, sim_xt, ...
        emulator_lambda,false);
    cov_xt_xt(:,:,ii) = cov_xt_xt(:,:,ii) + ...
        1e-5*eye(size(cov_xt_xt(:,:,ii)));
end

for kk = 1 : length(results)

    % Get posterior samples
    theta = results{kk}.theta1;

    % Get set up for loop that will find posterior predictive dist
    pred_xt_all = ([.5*ones(size(theta,1),size(xmin,2)) theta] - ...
        [zeros(1,size(xmin,2)) tmin])./[ones(1,size(xmin,2)) trange];
    prior_mean_xt_all = mean_sim(pred_xt_all) ; 

    % Get posterior predictive distribution for each output
    % Only do maxsimul samples at once, to ease computational expense
    maxsimul = 1000;
    for jj = 1 : ceil(size(theta,1)/maxsimul)
        startpt = (jj-1) * maxsimul + 1 ;
        endpt = min(jj * maxsimul, size(theta,1));
        pred_xt = pred_xt_all(startpt:endpt,:);
        for ii=1:2
            emulator_rho = results{kk}.settings.emulator_rho(:,ii);
            emulator_lambda = results{kk}.settings.emulator_lambda(ii);
            cov_pxt_xt = ...
                gp_cov(emulator_rho, pred_xt, sim_xt, ...
                emulator_lambda,false);
            cov_pxt_pxt = ...
                gp_cov(emulator_rho, pred_xt, pred_xt,...
                emulator_lambda,false);
            pred_mean_xt_current(:,ii) = ...
                prior_mean_xt_all(startpt:endpt,ii) + ...
                cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
                (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
            pred_cov_xt_current(:,:,ii) = ...
                cov_pxt_pxt - cov_pxt_xt * ...
                    (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
            fprintf('\n%g of %g, %g of 3 Done\n',jj,...
                ceil(size(theta,1)/maxsimul),ii);
            pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
            pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
                pred_cov_xt_current(:,:,ii);
        end
    end
    % Put back on original scale
    mean_y = results{kk}.settings.mean_y; 
    std_y = results{kk}.settings.std_y;
    pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
    mean(pred_mean_xt_os)


    % Now resample from the relevant gaussian dists, for full uncertainty
    for ii = 1:2
        pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
        pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
            pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
    end
    % Put back on original scale
    pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
    mean(pred_samps_xt_os)


    %%%%%%%%%%% Make a figure
    % xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
    % ylims = [0   95 ; 0   120 ; 0    3.5000];
    figure('Position',[10 10 600 300],'Color','white');
    for ii = 1:2
        subplot(1,3,ii);
%         histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
%             'EdgeColor','none');
%         hold on;
        histogram(pred_mean_xt_os(burn_in:end,ii),'Normalization','pdf',...
            'EdgeColor','none');
        xlabel(sprintf('y_%g',ii));
        set(gca,'Yticklabel',[]);
        set(gca,'Ytick',[]);
    %     xlim(xlims(ii,:));
    %     ylim(ylims(ii,:));
    drawnow;
    end

    % Save the posterior predictive output
    model_output.means = pred_mean_xt;
    model_output.vars = cell2mat(...
        arrayfun(@(b)diag(pred_cov_xt(:,:,b)),1:2,'UniformOutput',false));
    results{kk}.model_output = model_output;
end


% save(locstr,'results');

%% WTA Gather NSGA-II results for the new wind turbine applications
clc ; clearvars -except dpath ; close all ;

% Load the cost grid results
% locstr= ['C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\'...
%     'stored_data\2020-04-26_CTO_costgrid_size30'];
locstr = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\stored_data\2020-04-28_CTO_costgrid_size30';
load(locstr);

    
sim_t = results{1}.settings.sim_t1 .* results{1}.settings.range_t1 +...
    results{1}.settings.min_t1;

options = nsgaopt();                    % create default options structure
options.popsize = 50;                   % populaion size
options.maxGen  = 500;                  % max generation

options.numObj = 2;                     % number of objectives
options.numVar = size(results{1}.theta1,2); % number of design variables
options.numCons = 0;                    % number of constraints
options.lb = min(sim_t);               % lower bound of x
options.ub = max(sim_t);                % upper bound of x
options.objfun = @(x) emulator_mean(results{1},x);% objective function handle

options.plotInterval = 10;              % large interval for efficiency
options.outputInterval = 10;

tic;
nsga_result = nsga2(options);
elapsed_time = toc;
nsga_result.elapsed_time = elapsed_time;


% For convenience, add top-level array to result with final generation
% design settings and objective values.
final_theta = nan(options.popsize,options.numVar);
final_obj = nan(options.popsize,options.numObj);

for ii=1:options.popsize
    final_theta(ii,:) = nsga_result.pops(options.maxGen,ii).var;
    final_obj(ii,:) = nsga_result.pops(options.maxGen,ii).obj;
end

nsga_result.final_theta = final_theta;
nsga_result.final_obj = final_obj;

results{end+1}=nsga_result;
    


% Save NSGA-II result
% save(locstr,'results');

%% WTA CTO on bivariate objective over cost grid (for PF est) w/o obs error
clc ; clearvars -except dpath ; close all ;


%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
% Uncomment in order to reduce the size of the training set:
k=30; rng(355); raw_dat = raw_dat(sort(randsample(size(raw_dat,1),k)),:);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,[4,6]); 
clear raw_dat


% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Settings
obs_var = [1e-5 1e-9]; % Set observation error
obs_x_size  = 1; % Number of target observations
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 1e4; b = .5 ;

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;

% Emulator mean
%     mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 
% Evaluate with deg 2 poly regression
mod_terms = [];
coefs=[];
sim_x_01 = (sim_x - xmin)./xrange;
sim_t_01 = (sim_t - tmin)./trange;
sim_y_std = (sim_y - mean_y)./std_y;
for ii=1:size(sim_y,2)
    mod = polyfitn([sim_x_01 sim_t_01],sim_y_std(:,ii),2);
    mod_terms(:,:,ii) = mod.ModelTerms;
    coefs(:,ii) = mod.Coefficients;
end
mean_sim = @(xt_01) reshape(cell2mat(...
    arrayfun(@(row_idx) ...
    sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
    size(mod_terms,1),size(mod_terms,3)),1),(1:size(xt_01,1)),...
    'UniformOutput',false)'),size(xt_01,1),size(sim_y,2)) ;

% Define grid of targets
ngrid = 20 ; % Number of points in grid
% Minimum and maximum possible cost:
mincost = min(sim_y(:,2)); maxcost = max(sim_y(:,2)); 
des_obs_grid = [zeros(ngrid,1) linspace(96,352,ngrid)'] ;

% Perform CTO for each point in the grid
for ii = 1 : ngrid

    % Select the target design to be the point chosen in CTO (elsewhere)
    obs_y = repmat(des_obs_grid(ii,:),obs_x_size,1);

    % Get settings for DCTO
    settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
        obs_x,[],obs_y,...
        [],zeros(0,size(obs_y,2)),'min_x',xmin,'range_x',xrange,...
        'min_t1',tmin,'range_t1',trange,...
        'min_t2',[],'range_t2',[],...
        'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
        'obs_discrep',false,...
        'emulator_use',true,...
        'EmulatorMean',mean_sim,... 
        'EmulatorCovHypers','Default',...
        'des_discrep',false,...
        'obs_var',obs_var,...
        'obs_discrep_use_MLEs',false,...
        'doplot',doplot,...
        'verbose',verbose,...
        'obs_var_est',false,...
        'CTO',true);

    % Perform dual calibration
    % We need a loop because sometimes an ill-conditioned matrix early
    % in the burn-in makes the whole calibration fail.
    count = 0 ; err_count = 0 ; 
    while count == err_count
        try
            tic;
            res = MCMC_dual_calib(settings);
            elapsed_time = toc;
            res.elapsed_time = elapsed_time;
        catch ME
            warning('Warning: calibration failed. Retrying...');
            err_count = err_count + 1;
        end
        count = count + 1;
        if count >= 10, rethrow(ME) ; end
    end

    results{ii} = res;
    save temp

end


% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_CTO_costgrid_size',int2str(k)];
% save(locstr,'results');

%% WTA_alt CTO on 2 objectives over cost grid (PF est), no obs err var est.
clc ; clearvars -except dpath ; close all ;

% Store all results here:
all_results = {};

%%% Load raw data
locstr = [dpath,'stored_data\'...
    '2020-04-23_new_WTA_cases_size100'];
load(locstr);

for jj=1:length(dats)

    sim_x = dats{jj}.sim_x;
    sim_t = dats{jj}.sim_t;
    sim_y = dats{jj}.sim_y; 


    % Define inputs mins and ranges 
    xmin = min(sim_x);
    xrange = range(sim_x);
    tmin = min(sim_t);
    trange = range(sim_t);

    % Settings
    obs_x_size  = 1; % Number of target observations
    doplot = true;
    verbose = true;

    % Set number of draws, burn_in for each mcmc:
    M = 1e4; b = .5 ;

    % Get output means and stds
    mean_y = mean(sim_y) ; std_y = std(sim_y) ;

    % Get target design
    obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;

    % Emulator mean
%     mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 
    % Evaluate with deg 2 poly regression
    mod_terms = [];
    coefs=[];
    sim_x_01 = (sim_x - xmin)./xrange;
    sim_t_01 = (sim_t - tmin)./trange;
    sim_y_std = (sim_y - mean_y)./std_y;
    for ii=1:size(sim_y,2)
        mod = polyfitn([sim_x_01 sim_t_01],sim_y_std(:,ii),2);
        mod_terms(:,:,ii) = mod.ModelTerms;
        coefs(:,ii) = mod.Coefficients;
    end
    mean_sim = @(xt_01) reshape(cell2mat(...
        arrayfun(@(row_idx) ...
        sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
        size(mod_terms,1),size(mod_terms,3)),1),(1:size(xt_01,1)),...
        'UniformOutput',false)'),size(xt_01,1),size(sim_y,2)) ;

    % Define grid of targets
    ngrid = 20 ; % Number of points in grid
    % Minimum and maximum possible cost:
    mincost = min(sim_y(:,2)); maxcost = max(sim_y(:,2)); 
    mindefl = min(sim_y(:,1)); maxdefl = max(sim_y(:,1)); 
    des_obs_grid = [zeros(ngrid,1) linspace(mincost,maxcost,ngrid)'];
    obs_var_est_grid = repmat([false false],ngrid,1);
    obs_var_grid = repmat([1e-5 1e-9],ngrid,1);

    % Perform CTO for each point in the grid
    for ii = 1 : ngrid

        % Select the target design to be the point chosen in CTO (elsewhere)
        obs_y = repmat(des_obs_grid(ii,:),obs_x_size,1);

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
            obs_x,[],obs_y,...
            [],zeros(0,size(obs_y,2)),'min_x',xmin,'range_x',xrange,...
            'min_t1',tmin,'range_t1',trange,...
            'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',false,...
            'emulator_use',true,...
            'EmulatorMean',mean_sim,... 
            'EmulatorCovHypers','Default',...
            'des_discrep',false,...
            'obs_var',obs_var_grid(ii,:),...
            'obs_discrep_use_MLEs',false,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est_grid(ii,:),...
            'CTO',true);

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

        results{ii} = res;
        save temp

    end

    all_results{jj}=results;
    clear results
    
end

% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_3_CTO_grids_size',int2str(ngrid)];
% save(locstr,'all_results');

%% WTA_alt CTO, 2 objectives over 2-way grid (PF est), no obs err var est.
clc ; clearvars -except dpath ; close all ;

% Store all results here:
all_results = {};

%%% Load raw data
locstr = [dpath,'stored_data\'...
    '2020-04-23_new_WTA_cases_size100'];
load(locstr);

for jj=1:length(dats)

    sim_x = dats{jj}.sim_x;
    sim_t = dats{jj}.sim_t;
    sim_y = dats{jj}.sim_y; 


    % Define inputs mins and ranges 
    xmin = min(sim_x);
    xrange = range(sim_x);
    tmin = min(sim_t);
    trange = range(sim_t);

    % Settings
    obs_x_size  = 1; % Number of target observations
    doplot = false;
    verbose = false;

    % Set number of draws, burn_in for each mcmc:
    M = 1e4; b = .5 ;

    % Get output means and stds
    mean_y = mean(sim_y) ; std_y = std(sim_y) ;

    % Get target design
    obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;

    % Emulator mean
%     mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 
    % Evaluate with deg 2 poly regression
    mod_terms = [];
    coefs=[];
    sim_x_01 = (sim_x - xmin)./xrange;
    sim_t_01 = (sim_t - tmin)./trange;
    sim_y_std = (sim_y - mean_y)./std_y;
    for ii=1:size(sim_y,2)
        mod = polyfitn([sim_x_01 sim_t_01],sim_y_std(:,ii),2);
        mod_terms(:,:,ii) = mod.ModelTerms;
        coefs(:,ii) = mod.Coefficients;
    end
    mean_sim = @(xt_01) reshape(cell2mat(...
        arrayfun(@(row_idx) ...
        sum(coefs .* reshape(prod(xt_01(row_idx,:).^mod_terms,2),...
        size(mod_terms,1),size(mod_terms,3)),1),(1:size(xt_01,1)),...
        'UniformOutput',false)'),size(xt_01,1),size(sim_y,2)) ;

    % Define grid of targets
    ngrid = 20 ; % Number of points in grid
    % Minimum and maximum possible cost:
    mincost = min(sim_y(:,2)); maxcost = max(sim_y(:,2)); 
    mindefl = min(sim_y(:,1)); maxdefl = max(sim_y(:,1)); 
    des_obs_grid = [zeros(ngrid/2,1) linspace(mincost,maxcost,ngrid/2)'; ...
        linspace(mindefl,maxdefl,ngrid/2)' zeros(ngrid/2,1)] ;
    obs_var_est_grid = [repmat([false false],ngrid/2,1);...
        repmat([false false],ngrid/2,1)];
    obs_var_grid = [repmat([1e-5 1e-9],ngrid/2,1);...
        repmat([1e-9 1e-5],ngrid/2,1)]; % initial guess of obs error var

    % Perform CTO for each point in the grid
    for ii = 1 : ngrid

        % Select the target design to be the point chosen in CTO (elsewhere)
        obs_y = repmat(des_obs_grid(ii,:),obs_x_size,1);

        % Get settings for DCTO
        settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
            obs_x,[],obs_y,...
            [],zeros(0,size(obs_y,2)),'min_x',xmin,'range_x',xrange,...
            'min_t1',tmin,'range_t1',trange,...
            'min_t2',[],'range_t2',[],...
            'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
            'obs_discrep',false,...
            'emulator_use',true,...
            'EmulatorMean',mean_sim,... 
            'EmulatorCovHypers','Default',...
            'des_discrep',false,...
            'obs_var',obs_var_grid(ii,:),...
            'obs_discrep_use_MLEs',false,...
            'doplot',doplot,...
            'verbose',verbose,...
            'obs_var_est',obs_var_est_grid(ii,:),...
            'CTO',true);

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

        results{ii} = res;
        save temp

    end

    all_results{jj}=results;
    clear results
    
end

% locstr = [dpath,'stored_data\'...
%     '2019-11-05_CTO_costgrid'];
locstr = [dpath,'stored_data\'...
    datestr(now,'yyyy-mm-dd'),'_3_CTO_double_grid_size',int2str(ngrid)];
save(locstr,'all_results');