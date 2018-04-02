% Figures for fourth exam


clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');

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


%% Comparison of different versions of observation variance
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
    
    subplot(1,3,2)
    plot(bsamps(:,2),'ko');
    title('Thickness draws');
    xlabel('Draw'); ylabel('Thickness');
    
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
    plot(rr(:,1),rr(:,2),'r.');
    hold off;
    title({'(Logit) VF vs. thickness', 'with proposal covariance'});
    xlabel('Logit VF'); ylabel('Logit thickness');
    suptitle(strcat('Desired data: deflection',{' '},...
        num2str(desired_data(jj,1)),{', '},...
        'cost', {' '}, num2str(desired_data(jj,2))));
    
    scale = 0.1;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
    
    waitforbuttonpress;
    %saveas(h,...
    %    sprintf('FIG%d.png',jj));
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

h = figure('rend','painters','pos',[10 10 1000 250]);
subplot(1,2,1);
acf(bsamps1(:,1),50);
ylim([-0.1 1]);
title('Boundary conditions present');
scale = 0.15;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)
subplot(1,2,2);
acf(bsamps2(:,1),50);
ylim([-0.1 1]);
title('Boundary conditions removed')
scale = 0.15;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)
suptitle('Volume fraction autocorrelation');
saveas(h,sprintf('ACF_bnd_cnds_fig.png'));

