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

load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\'...
    '30_MCMCs2.mat']);
results = cell(2,1);
results{1}=all_results2{2}{2};
results{2}=all_results2{1}{3};

for jj = 1:2
    bsamps = results{jj}.samples(2002:end,:);
    bsigs = results{jj}.sigma2(2002:end,:);
    
    subplot(1,3,1)
    plot(bsamps(:,1),'ko');
    
    subplot(1,3,2)
    plot(bsamps(:,2),'ko');
    
    subplot(1,3,3)
    plot(logit(bsamps(:,1)),logit(bsamps(:,2)),'ko');
    mu = mean(logit(bsamps));
    Sig = results{jj}.Sigma;
    rr = mvnrnd(mu,Sig,250);
    hold on;
    plot(rr(:,1),rr(:,2),'r.');
    hold off;
    
    waitforbuttonpress;
end


