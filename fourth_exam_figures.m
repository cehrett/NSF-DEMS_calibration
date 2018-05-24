% Figures for fourth exam




%% Clear and add paths
clc; clear all; close all;

addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');

% Set path string
direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

%% Comparison of different versions of observation variance
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    'results_1d_homosked6sdObs.mat']);
results1=results;
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    'results_1d_homosked2sdposObs.mat']);
results2=results;
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    'results_1d_heteroskedPriorObs.mat']);
results3=results;
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    'results_1d_homoskedPriorObs.mat']);
results4=results;


bsamps1 = results1.samples(2000:end,:);
bsamps2 = results2.samples(2000:end,:);
bsamps4 = results3.samples(2000:end,:);
bsamps3 = results4.samples(2000:end,:);

subplot(2,2,1);
plot(bsamps1(:,1),bsamps1(:,2),'ko');
subplot(2,2,2);
plot(bsamps2(:,1),bsamps2(:,2),'ko');
subplot(2,2,3);
plot(bsamps3(:,1),bsamps3(:,2),'ko');
subplot(2,2,4);
plot(bsamps4(:,1),bsamps4(:,2),'ko');

subplot(2,2,1);
plot(bsamps1(:,1),'ko');
title({'Homoskedastic version,','constant \sigma^2 = 6'});
xlim([0,8000])
subplot(2,2,2);
plot(bsamps2(:,1),'ko');
title({'Heteroskedastic version,','constant 2 s.d.''s positive'});
xlim([0,8000])
subplot(2,2,3);
plot(bsamps3(:,1),'ko');
title({'Homoskedastic version,','1/\sigma^2 prior'});
xlim([0,8000])
subplot(2,2,4);
plot(bsamps4(:,1),'ko');
title({'Heteroskedastic version,','1/\sigma^2 prior'});
xlim([0,8000])


%% Comparison of 1d vs 0d
clc; clear all; close all;

load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    'results_1d_heteroskedPriorObs.mat']);
results1=results;
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    'results_0d_heteroskedPriorObs.mat']);
results2=results;

bsamps1 = results1.samples(2000:end,:);
bsamps2 = results2.samples(2000:end,:);
bsigs1 = results1.sigma2(2000:end,:);
bsigs2 = results2.sigma2(2000:end,:);

subplot(2,5,1);
plot(bsamps1(:,1),'ko');
xlim([0,8000])
subplot(2,5,2);
plot(bsamps1(:,2),'ko');
xlim([0,8000])
subplot(2,5,3);
plot(bsigs1(:,1),'ko');
xlim([0,8000]);
subplot(2,5,4);
plot(bsigs1(:,2),'ko');
xlim([0,8000]);
subplot(2,5,5);
plot(bsigs1(:,3),'ko');
xlim([0,8000]);
subplot(2,5,6)
plot(bsamps2(:,1),'ko');
xlim([0,8000])
subplot(2,5,7);
plot(bsamps2(:,2),'ko');
xlim([0,8000])
subplot(2,5,8);
plot(bsigs2(:,1),'ko');
xlim([0,8000]);
subplot(2,5,9);
plot(bsigs2(:,2),'ko');
xlim([0,8000]);
subplot(2,5,10);
plot(bsigs2(:,3),'ko');
xlim([0,8000]);

%% Example of healthy and unhealthy results from current version
clc; clear all; close all;
logit = @(x) log(x./(1-x));

load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
    '30_MCMCs.mat']);
results = cell(2,1);
results{1}=all_results2{2}{2};
results{2}=all_results2{1}{3};
desired_data = [ results{1}.desired_data ; results{2}.desired_data ];

h = figure('rend','painters','pos',[10 10 1000 275]);
for jj = 1:2
    bsamps = results{jj}.samples(2002:end,:);
    bsigs = results{jj}.sigma2(2002:end,:);
    
    subplot(1,3,1)
    plot(bsamps(:,1),'ko');
    title('Volume fraction draws');
    xlabel('Draw'); ylabel('VF');
    
    scale = 0.1;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
    ylim([0,1]);
    
    subplot(1,3,2)
    plot(bsamps(:,2),'ko');
    title('Thickness draws');
    xlabel('Draw'); ylabel('Thickness');
    ylim([0,1]);
    
    scale = 0.1;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
    
    subplot(1,3,3)
    plot(logit(bsamps(:,1)),logit(bsamps(:,2)),'ko');
    mu = mean(logit(bsamps));
    Sig = results{jj}.Sigma;
    rr = mvnrnd(mu,Sig,250);
    hold on;
    %plot(rr(:,1),rr(:,2),'r.');
    hold off;
    title({'(Logit) VF vs. thickness'});
    xlabel('Logit VF'); ylabel('Logit thickness');
    suptitle(strcat('Desired data: deflection',{' '},...
        num2str(round(desired_data(jj,1),2)),{', '},...
        'cost', {' '}, num2str(round(desired_data(jj,2)))));
    
    scale = 0.1;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
    
    waitforbuttonpress;
    saveas(h,...
        sprintf('FIG%d_posterior_samps_at_spec_cost.png',jj));
end


%% Comparison w/ and w/o boundary constraints

clc; clear all; close all;

load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
   'results_0d_heteroskedPriorObs_hasBndCnds.mat']);
results1 = results;
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
   'results_0d_heteroskedPriorObs.mat']);
results2 = results;

bsamps1 = results1.samples(2002:end,:);
bsamps2 = results2.samples(2002:end,:);

h = figure('rend','painters','pos',[10 10 520 250]);
%subplot(1,2,1);
bs1acf=acf(bsamps1(:,1),50);

% title('Boundary conditions present');
% scale = 0.15;
% pos = get(gca, 'Position');
% pos(2) = pos(2)+scale*pos(4);
% pos(4) = (1-scale)*pos(4);
% set(gca, 'Position', pos)
% subplot(1,2,2);
bs2acf=acf(bsamps2(:,1),50);
% ylim([-0.1 1]);
% title('Boundary conditions removed')
% scale = 0.15;
% pos = get(gca, 'Position');
% pos(2) = pos(2)+scale*pos(4);
% pos(4) = (1-scale)*pos(4);
% set(gca, 'Position', pos)
plot(bs1acf,'-o');
hold on;
plot(bs2acf,'-o');
ylim([-0.05 1]);
xlim([1,50]);
line([1 50] , [ 0 0 ],'Color','k','LineWidth',1.5);
xlabel('Lag'); ylabel('Autocorrelation');
title('Volume fraction autocorrelation');
legend('Without boundary constraints','With boundary constraints');
saveas(h,sprintf('ACF_bnd_cnds_fig.png'));

%% Lambda_cost pareto bands (original)
load([dpath,'stored_data\'...
    'results_Cost_lambda_grid_exploration'],...
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
    pdo(ii,:) = quantile(results{ii}.model_output.by_sample,0.5);
    pso(ii,:) = norminv(1-alpha/2) * mean(results{ii}.model_output.sds);
    plo(ii,:) = quantile(results{ii}.model_output.by_sample,alpha/2);
    puo(ii,:) = quantile(results{ii}.model_output.by_sample,1-alpha/2);
    cost_lambda(ii) = results{ii}.Cost_lambda;
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

% Begin the figure
h=figure('rend','painters','pos',[10 10 1200 400]);
x = 0:1:100;
subplot(1,3,1)
% Get main curve
pcost = pchip(cost_lambda,post_cost_median,x);
% Get upper and lower 0.05 quantiles curves
pcostuq = pchip(cost_lambda,post_cost_uq,x);
pcostlq = pchip(cost_lambda,post_cost_lq,x);
f=fill([ x , fliplr(x) ], [pcostuq, fliplr(pcostlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(...%cost_lambda,post_cost_median,'r',...
    x,pcost,'-r',...
    x,pcostuq,'-k',...
    x,pcostlq,'-k');
% Plot 2sd errbar
% errorbar(cost_lambda,post_cost_mean,post_cost_sd,'ob'); 
xl1=xlabel('\lambda_c_o_s_t');
ylim(ylim_cost);
ylabel('Cost');


subplot(1,3,2)
% Get main curve
pdefl = pchip(cost_lambda,post_defl_median,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(cost_lambda,post_defl_uq,x);
pdefllq = pchip(cost_lambda,post_defl_lq,x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost_lambda,post_defl_median,'r',...
    x,pdefl,'-r',...
    x,pdefluq,'-k',...
    x,pdefllq,'-k');
xl2=xlabel('\lambda_c_o_s_t');
ylim(ylim_defl);
ylabel('Deflection');

subplot(1,3,3)
% Get main curve
protn = pchip(cost_lambda,post_rotn_median,x);
% Get upper and lower 0.05 quantiles curves
protnuq = pchip(cost_lambda,post_rotn_uq,x);
protnlq = pchip(cost_lambda,post_rotn_lq,x);
f=fill([ x , fliplr(x) ], [protnuq, fliplr(protnlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(cost_lambda,post_rotn_median,'r',...
    x,protn,'-r',...
    x,protnuq,'-k',...
    x,protnlq,'-k');
xl3=xlabel('\lambda_c_o_s_t');
ylim(ylim_rotn);
ylabel('Rotation');


% Save the figure temporarily so we can mess with it later, because
% suptitle seems to mess things up somehow for making changes after calling
% it
savefig(h,'tempfig');

% Put a main title over anything, and fix any misplaced labels etc
suptitle(['Performance metrics vs. \lambda_c_o_s_t,',...
    ' with ',num2str(cred_level),'% credible interval']); 
p = get(xl1,'position');
set(xl1,'position',p + [0 6 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.003 0])
p = get(xl3,'position');
set(xl3,'position',p + [0 0.00045 0])
figpos = get(h,'pos');

% Save the figure
saveas(h,'FIG_cost_lambda.png');

% Now add in code uncertainty. That is, the above assumes that the GP
% emulator nails the FE code precisely. But of course the GP emulator has
% nonnegligible variance. That's the code uncertainty. So our confidence
% bands should reflect it. So we add it in here, by dropping the
% appropriate multiple of the sd from each lower quantile and adding it to
% each upper quantile.
% First, open the figure prior to calling suptitle.
h=openfig('tempfig');
subplot(1,3,1);
pcostuq_code_uncert = pchip(cost_lambda,post_cost_uq_cu,x);
pcostlq_code_uncert = pchip(cost_lambda,post_cost_lq_cu,x);
f=fill([ x , fliplr(x) ], [pcostuq_code_uncert,...
    fliplr(pcostuq)],'b');
ff=fill([ x , fliplr(x) ], [pcostlq_code_uncert,...
    fliplr(pcostlq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
ylim(ylim_cost);
xl1=xlabel('\lambda_c_o_s_t');

subplot(1,3,2);
pdefluq_code_uncert = pchip(cost_lambda,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(cost_lambda,post_defl_lq_cu,x);
f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
    fliplr(pdefluq)],'b');
ff=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
    fliplr(pdefllq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
ylim(ylim_defl);
xl2=xlabel('\lambda_c_o_s_t');

subplot(1,3,3);
protnuq_code_uncert = pchip(cost_lambda,post_rotn_uq_cu,x);
protnlq_code_uncert = pchip(cost_lambda,post_rotn_lq_cu,x);
f=fill([ x , fliplr(x) ], [protnuq_code_uncert,...
    fliplr(protnuq)],'b');
ff=fill([ x , fliplr(x) ], [protnlq_code_uncert,...
    fliplr(protnlq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
ylim(ylim_rotn);
xl3=xlabel('\lambda_c_o_s_t');

suptitle(['Performance metrics vs. \lambda_c_o_s_t,',...
    ' with ',num2str(cred_level),'% credible interval '...
    'including code uncertainty']); 
set(h,'pos',figpos);
p = get(xl1,'position');
set(xl1,'position',p + [0 6 0]);
p = get(xl2,'position');
set(xl2,'position',p + [0 0.003 0])
p = get(xl3,'position');
set(xl3,'position',p + [0 0.00045 0])

saveas(h,'FIG_cost_lambda_code_uncert.png');

delete('tempfig.fig');

%% Lambda_cost pareto bands (updated)
% Start out by getting the plot of Cost vs lambda_cost
load([dpath,'stored_data\'...
    'results_Cost_lambda_grid_exploration'],...
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
    pdo(ii,:) = quantile(results{ii}.model_output.by_sample,0.5);
    pso(ii,:) = norminv(1-alpha/2) * mean(results{ii}.model_output.sds);
    plo(ii,:) = quantile(results{ii}.model_output.by_sample,alpha/2);
    puo(ii,:) = quantile(results{ii}.model_output.by_sample,1-alpha/2);
    cost_lambda(ii) = results{ii}.Cost_lambda;
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

% Begin the figure
h=figure('rend','painters','pos',[10 10 1200 400]);
x = 0:1:100;
subplot(1,3,1)
% Get main curve
pcost = pchip(cost_lambda,post_cost_median,x);
% Get upper and lower 0.05 quantiles curves
pcostuq = pchip(cost_lambda,post_cost_uq,x);
pcostlq = pchip(cost_lambda,post_cost_lq,x);
f=fill([ x , fliplr(x) ], [pcostuq, fliplr(pcostlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(...%cost_lambda,post_cost_median,'r',...
    x,pcost,'-r',...
    x,pcostuq,'-k',...
    x,pcostlq,'-k');
xl1=xlabel('\lambda_c_o_s_t');
ylim(ylim_cost);
ylabel('Cost');
load([dpath,'stored_data\'...
    'results_Cost_lambda_grid_exploration'],...
    'results');

% Collect result outputs in one place
outputs=[];
sds    =[];
for ii = 1:size(results,1)
    outputs = [outputs ; results{ii}.model_output.by_sample ] ;
    sds     = [ sds    ; results{ii}.model_output.sds       ] ;
end

mesh=10; % This describes how fine we want the cost grid to 
meshin=1/mesh; % Useful
lowval = round(min(outputs(:,3))/mesh)*mesh; % round to multiple of mesh
highval = round(max(outputs(:,3))/mesh)*mesh;
vals = lowval:mesh:highval;

% Now separate the outputs by rounding to nearest point in cost grid
grid_outputs = cell(length(vals),1);
grid_sds = cell(length(vals),1);
for ii = 1 : length(vals)
    idx = (round(outputs(:,3)/mesh)*mesh == vals(ii));
    grid_outputs{ii} = outputs(idx,:);
    grid_sds{ii}     = sds(idx,:)    ;
end
    
% Collect costs, and posterior mean and quantiles for costs, defl, rot
m=size(grid_outputs,1); % Store number of target cost_lambdas
cred_level = 90; % Set desired level for credible bands (in %)
alpha = (100-cred_level)/100; % Convert cred_level to alpha level
pmo = zeros(m,3); % This will store posterior mean output of emulator
pdo = zeros(m,3); % ``'' median output
pso = zeros(m,3); % ``'' appropriate multiple of standard deviations
plo = zeros(m,3); % ``'' lower (alpha/2) quantile
puo = zeros(m,3); % ``'' upper (alpha/2) quantile
for ii = 1:m % This loop populates the above arrays
    pdo(ii,:) = quantile(grid_outputs{ii},0.5);
    pso(ii,:) = norminv(1-alpha/2) * mean(grid_sds{ii});
    plo(ii,:) = quantile(grid_outputs{ii},alpha/2);
    puo(ii,:) = quantile(grid_outputs{ii},1-alpha/2);
end
% Now we break the arrays up each into 3 vectors, one for each output
%post_cost_median = pdo(:,3);
post_defl_median = pdo(:,1);
post_rotn_median = pdo(:,2);
%post_cost_sd = pso(:,3);
post_defl_sd = pso(:,1);
post_rotn_sd = pso(:,2);
%post_cost_lq = plo(:,3);
%post_cost_uq = puo(:,3);
post_defl_lq = plo(:,1);
post_defl_uq = puo(:,1);
post_rotn_lq = plo(:,2);
post_rotn_uq = puo(:,2);
% Get quantiles plus code uncertainty
%post_cost_uq_cu = post_cost_uq + post_cost_sd;
%post_cost_lq_cu = post_cost_lq - post_cost_sd;
post_defl_uq_cu = post_defl_uq + post_defl_sd;
post_defl_lq_cu = post_defl_lq - post_defl_sd;
post_rotn_uq_cu = post_rotn_uq + post_rotn_sd;
post_rotn_lq_cu = post_rotn_lq - post_rotn_sd;
% Get ylims for the two sets of plots
ylimrat=1.01;
%ylim_cost = [min(post_cost_lq_cu)/ylimrat max(post_cost_uq_cu)*ylimrat];
ylim_defl = [min(post_defl_lq_cu)/ylimrat max(post_defl_uq_cu)*ylimrat];
ylim_rotn = [min(post_rotn_lq_cu)/ylimrat max(post_rotn_uq_cu)*ylimrat];


% Begin the figure
%h=figure('rend','painters','pos',[10 10 1200 400]);

% subplot(1,3,1)
% % Get main curve
% pcost = pchip(cost_lambda,post_cost_median,x);
% % Get upper and lower 0.05 quantiles curves
% pcostuq = pchip(cost_lambda,post_cost_uq,x);
% pcostlq = pchip(cost_lambda,post_cost_lq,x);
% f=fill([ x , fliplr(x) ], [pcostuq, fliplr(pcostlq)],'k');
% set(f,'facealpha',.25);
% hold on;
% plot(cost_lambda,post_cost_median,'or',...
%     x,pcost,'-r',...
%     x,pcostuq,'-k',...
%     x,pcostlq,'-k');
% % Plot 2sd errbar
% % errorbar(cost_lambda,post_cost_mean,post_cost_sd,'ob'); 
% xl1=xlabel('\lambda_c_o_s_t');
% ylim(ylim_cost);
% ylabel('Cost');

% Now make the two new plots: defl vs cost, rotn vs cost:
x = min(vals):1:max(vals);
subplot(1,3,2)
% Get main curve
pdefl = pchip(vals,post_defl_median,x);
% Get upper and lower 0.05 quantiles curves
pdefluq = pchip(vals,post_defl_uq,x);
pdefllq = pchip(vals,post_defl_lq,x);
f=fill([ x , fliplr(x) ], [pdefluq, fliplr(pdefllq)],'k');
set(f,'facealpha',.25);
hold on;
plot(vals,post_defl_median,'r',...
    x,pdefl,'-r',...
    x,pdefluq,'-k',...
    x,pdefllq,'-k');
xl2=xlabel('Cost ($)');
ylim(ylim_defl);
ylabel('Deflection');

subplot(1,3,3)
% Get main curve
protn = pchip(vals,post_rotn_median,x);
% Get upper and lower 0.05 quantiles curves
protnuq = pchip(vals,post_rotn_uq,x);
protnlq = pchip(vals,post_rotn_lq,x);
f=fill([ x , fliplr(x) ], [protnuq, fliplr(protnlq)],'k');
set(f,'facealpha',.25);
hold on;
plot(vals,post_rotn_median,'r',...
    x,protn,'-r',...
    x,protnuq,'-k',...
    x,protnlq,'-k');
xl3=xlabel('Cost ($)');
ylim(ylim_rotn);
ylabel('Rotation');

% Note to get plot one you have to use the previous figure section above.


% Save the figure temporarily so we can mess with it later, because
% suptitle seems to mess things up somehow for making changes after calling
% it
figpos = get(h,'pos');
savefig(h,'tempfig');

% Put a main title over anything, and fix any misplaced labels etc
suptitle(['Cost vs. \lambda_c_o_s_t, and performance metric vs. cost',...
    ' with ',num2str(cred_level),'% credible interval']); 
p = get(xl1,'position');
set(xl1,'position',p + [0 4 0]);

% Save the figure
saveas(h,'FIG_cost_lambda_upd.png');

% Now add in code uncertainty. That is, the above assumes that the GP
% emulator nails the FE code precisely. But of course the GP emulator has
% nonnegligible variance. That's the code uncertainty. So our confidence
% bands should reflect it. So we add it in here, by dropping the
% appropriate multiple of the sd from each lower quantile and adding it to
% each upper quantile.
% First, open the figure prior to calling suptitle.
h=openfig('tempfig');
subplot(1,3,1);
x = 0:1:100;
pcostuq_code_uncert = pchip(cost_lambda,post_cost_uq_cu,x);
pcostlq_code_uncert = pchip(cost_lambda,post_cost_lq_cu,x);
f=fill([ x , fliplr(x) ], [pcostuq_code_uncert,...
    fliplr(pcostuq)],'b');
ff=fill([ x , fliplr(x) ], [pcostlq_code_uncert,...
    fliplr(pcostlq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
ylim(ylim_cost);
xl1=xlabel('\lambda_c_o_s_t');

subplot(1,3,2);
x = min(vals):1:max(vals);
pdefluq_code_uncert = pchip(vals,post_defl_uq_cu,x);
pdefllq_code_uncert = pchip(vals,post_defl_lq_cu,x);
f=fill([ x , fliplr(x) ], [pdefluq_code_uncert,...
    fliplr(pdefluq)],'b');
ff=fill([ x , fliplr(x) ], [pdefllq_code_uncert,...
    fliplr(pdefllq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
ylim(ylim_defl);
xl2=xlabel('Cost ($)');

subplot(1,3,3);
protnuq_code_uncert = pchip(vals,post_rotn_uq_cu,x);
protnlq_code_uncert = pchip(vals,post_rotn_lq_cu,x);
f=fill([ x , fliplr(x) ], [protnuq_code_uncert,...
    fliplr(protnuq)],'b');
ff=fill([ x , fliplr(x) ], [protnlq_code_uncert,...
    fliplr(protnlq)],'b');
set(f,'facealpha',.25);
set(ff,'facealpha',.25);
ylim(ylim_rotn);
xl3=xlabel('Cost ($)');

suptitle(['Cost vs. \lambda_c_o_s_t, and performance metrics vs. cost',...
    ' with ',num2str(cred_level),'% credible interval '...
    'including code uncertainty']); 
set(h,'pos',figpos);
p = get(xl1,'position');
set(xl1,'position',p + [0 4 0]);
% p = get(xl2,'position');
% set(xl2,'position',p + [0 0.003 0])
% p = get(xl3,'position');
% set(xl3,'position',p + [0 0.00045 0])

saveas(h,'FIG_cost_lambda_code_uncert_upd.png');

delete('tempfig.fig');

%% Posterior theta distribution over lambda_cost grid
% Load data

% Let's scatterplotclc; clearvars -except dpath ; close all;

load([dpath,'stored_data\'...
    'results_Cost_lambda_grid_exploration'],...
    'results');

% Get data in one big array, [theta lambda_cost]
dat = [];
for ii = 1:size(results,1)
    sos = results{ii}.samples_os(results{ii}.burn_in:end,:);
    dat = [dat ; sos ...
        results{ii}.Cost_lambda * ones(size(sos,1),1) ];
end

colors = [dat(dat(:,3)<10,3) * 5 ; dat(dat(:,3)>10,3)+41];
scatter3(dat(:,1),dat(:,2),dat(:,3),5,colors); axis vis3d;

% Now make a scatterplot only of lower lambda_cost values. Rather than set
% the height proportional to lambda_cost, just set them at regular
% intervals and label the ticks with the real values. To do this, we need
% to make a vector of the heights at regular intervals.
upto=25;
scatter3(dat(dat(:,3)<upto,1),dat(dat(:,3)<upto,2),dat(dat(:,3)<upto,3)...
    ,5,colors(dat(:,3)<upto)); axis vis3d;