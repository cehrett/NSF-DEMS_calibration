% Dual calibration figures
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

%% Show example objective function
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
f=figure();
xidx=100;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
surf(tt1 * t1range + t1min,tt2*t2range+t2min,reshape(Y(:,xidx,:),100,100));

% Get a rotating gif
set(f,'color','white');
title('Objective function at x=1');
xlabel('\theta_1'); ylabel('\theta_2'); zlabel('f(x,\theta_1,\theta_2)');
viewpt = [-27.2667 10.4000];
view(viewpt);
axis vis3d;
% gif('FIG_dual_calib_obf_fn.gif','frame',gcf);
nfms = 120;
for ii = 1:nfms
    viewpt = viewpt + [ 360/nfms 0 ];
    view(viewpt);
    pause(.01);
%     gif
end

% Now get a view of it 2D, f vs theta1
a1 = gca;
f2 = figure();
copyobj(a1,f2);
set(f2,'color','white');
view([0,0]); hold on;
plot3([2 2],[0 0],[0 1],'r','LineWidth',2);
label = sprintf('True \\theta_1\nvalue');
text(2.05,0,0.2,label);
% saveas(f2,'FIG_dual_calib_true_theta1.png');

% Now get another 2D view, this time f vs theta2 at true theta1 value
f3 = figure('pos',[10 10 300 250]);
theta2 = tt2 * t2range + t2min;
theta2 = mean(theta2(17:18,:));
y = reshape(Y(:,xidx,:),100,100);
y = mean(y(17:18,:));
% Save the data being plotted for use in other figures
ytheta = struct('theta2',theta2,'y',y);
% save([dpath,'dual_calib\dual_calib_stored_data\'...
%     '2019-01-30_dual_calib_theta2_vs_y'],...
%     'ytheta');
plot(theta2,y,'LineWidth',2);
set(gca, 'XTick', [0 4/3 5]);
hold on;
plot([4/3 4/3],get(gca,'YLim'),'r','LineWidth',2);
label = sprintf('Optimal\n\\theta_2 value');
text(4/3+.1,.9,label);
title('Objective function at x=1, \theta_1=2');
xlabel('\theta_2'); ylabel('f(x,\theta_1,\theta_2)');
set(f3,'color','white');
% saveas(f3,'FIG_dual_calib_optimal_theta2.png');

%% Show posterior distributions of calibration parameters with priors
clc ; clearvars -except dpath ; close all ;

% Load results
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-02-08_dual_calib_nodiscr']);
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
    'EdgeColor','none','FaceColor','b','FaceAlpha',.9,'BinWidth',0.03);
% Plot true theta1
plot([2 2],get(gca,'YLim'),'--r','LineWidth',1.5);
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
    'EdgeColor','none','FaceColor','b','FaceAlpha',.9);
% Plot true theta2
yyaxis left ;
plot([4/3 4/3],get(gca,'YLim'),'--r','LineWidth',1.5);
% Put a legend on it
lg2 = legend('Prior dist.','Posterior dist.','Optimal value');
title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
% Save it as is:
set(f2,'color','white');
% saveas(f2,'FIG_dual_calib_post_theta2-2.png');

% Now add to the most recent plot:
% Add the shape of the objective function at true value of theta1
yyaxis right;
ylim([1 1.6]);
load([dpath,'dual_calib\dual_calib_stored_data\'...
    '2019-01-30_dual_calib_theta2_vs_y'],...
    'ytheta');
plot(ytheta.theta2,1./ytheta.y,'Color',right_color,'LineWidth',1.5);
lg2.String{4} = 'Inverse obj. fn.';
% Save it:
% saveas(f2,'FIG_dual_calib_post_theta2_w_obj_fn.png');