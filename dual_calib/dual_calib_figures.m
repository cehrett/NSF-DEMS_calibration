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

%% Examine four different discrepancy versions
clc ; clearvars -except dpath ; close all ;

f=figure('pos',[20 20 700 600]);
set(f,'color','white');

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

%%% Loop through all discrepancies and plot each
for ii=1:4
    subplot(2,2,ii);
    discrep = ii-1 ; % Select which discrepancy
    Y = reshape(...
        dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
        T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
    Yd= reshape(...
        dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
        T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

    % Take a look
    xidx=50;
    xx=reshape(X(:,xidx,:),100,100);
    % Get value of x
    xval = xx(1)*xrange + xmin;
    tt1=reshape(T1(:,xidx,:),100,100);
    tt2=reshape(T2(:,xidx,:),100,100);
    Discrep = Yd-Y;
    ea=.5;
    surf(tt1*t1range+t1min,tt2*t2range+t2min,...
        reshape(Y(:,xidx,:),100,100),...
        'EdgeAlpha',ea);
    hold on;
    surf(tt1*t1range+t1min,tt2*t2range+t2min,...
        reshape(Yd(:,xidx,:),100,100),...
        'EdgeAlpha',ea);
    surf(tt1*t1range+t1min,tt2*t2range+t2min,...
        reshape(Discrep(:,xidx,:),100,100),...
        'EdgeAlpha',ea);

    
%     title(sprintf('Objective function at x=%0.2f',xval));
    xlabel('\theta_1');ylabel('\theta_2');zlabel('f(x,\theta_1,\theta_2)');
    axis vis3d;
    view([-110.4000    6.5334]);
end

%%% Fix sizing
f.Children(4).Position = [0.0 0.505 0.575 0.55];
f.Children(3).Position = [0.5 0.505 0.575 0.55];
f.Children(2).Position = [0.0 0.005 0.575 0.55];
f.Children(1).Position = [0.5 0.005 0.575 0.55];

% % Get a rotating gif
% viewpt = [-27.2667 10.4000];
% view(viewpt);
% % gif('FIG_dual_calib_obf_fn.gif','frame',gcf);
% nfms = 120;
% for ii = 1:nfms
%     viewpt = viewpt + [ 360/nfms 0 ];
%     view(viewpt);
%     pause(.01);
% %     gif
% end

% Save it:
% saveas(f,'FIG_four_discrepancies.png');

%% Examine six different discrepancy versions
clc ; clearvars -except dpath ; close all ;

f=figure('pos',[640 5 540 800]);
set(f,'color','white');

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

discrep_title_content = [ {1;'c = 1.5'}, {1;'c = 3.5'}, ...
    {2;'c = .15, d = 0.075'}, {2;'c = .65, d = 0.075'}, ...
    {3;'c = 0.055, d=0'}, {3;'c = 0.055, d = 0.1'} ] ;

%%% Loop through all discrepancies and plot each
for ii=1:6
    subplot(3,2,ii);
    discrep = ii ; % Select which discrepancy
    Y = reshape(...
        dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
        T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
    Yd= reshape(...
        dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
        T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

    % Take a look
    xidx=50;
    xx=reshape(X(:,xidx,:),100,100);
    % Get value of x
    xval = xx(1)*xrange + xmin;
    tt1=reshape(T1(:,xidx,:),100,100);
    tt2=reshape(T2(:,xidx,:),100,100);
    Discrep = Yd-Y;
    ea=.25;
    surf(tt1*t1range+t1min,tt2*t2range+t2min,...
        reshape(Y(:,xidx,:),100,100),...
        'EdgeAlpha',ea);
    hold on;
    surf(tt1*t1range+t1min,tt2*t2range+t2min,...
        reshape(Yd(:,xidx,:),100,100),...
        'EdgeAlpha',ea);
    surf(tt1*t1range+t1min,tt2*t2range+t2min,...
        reshape(Discrep(:,xidx,:),100,100),...
        'EdgeAlpha',ea);
    
    % Sort out title and labels
    dtc = discrep_title_content(:,ii);
    title(sprintf('Discrepancy %d, %s',dtc{:}));
    xlabel('t_1');ylabel('t_2');zlabel('f(x,t_1,t_2)');
    axis vis3d;
    view([-110.4000    6.5334]);
end

%%% Fix sizing
f.Children(6).Position = [0.0 0.685 0.575 0.32];
f.Children(5).Position = [0.5 0.685 0.575 0.32];
f.Children(4).Position = [0.0 0.355 0.575 0.32];
f.Children(3).Position = [0.5 0.355 0.575 0.32];
f.Children(2).Position = [0.0 0.025 0.575 0.32];
f.Children(1).Position = [0.5 0.025 0.575 0.32];

% % Get a rotating gif
% viewpt = [-27.2667 10.4000];
% view(viewpt);
% % gif('FIG_dual_calib_obf_fn.gif','frame',gcf);
% nfms = 120;
% for ii = 1:nfms
%     viewpt = viewpt + [ 360/nfms 0 ];
%     view(viewpt);
%     pause(.01);
% %     gif
% end

% Save it:
% saveas(f,'FIG_six_discrepancies.png');

%% Results from calibration with various discrepancies
clc ; clearvars -except dpath ; close all ;

% Load results
discrep = 6; % Change this to whichever discrepancy is desired
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-15_dual_calib_discrep',int2str(discrep)]);
load(locstr,'results');


% Define true theta1
theta1 = 2;

burn_in = results.settings.burn_in;
min_t1 = results.settings.min_t1; range_t1 = results.settings.range_t1;
min_t2 = results.settings.min_t2; range_t2 = results.settings.range_t2;

% First, get prior and posterior theta1
f = figure('pos',[10 10 600 250]);
subplot(1,2,1);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 with true value marked
histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',.85,'BinWidth',0.075);
% Plot true theta1
set(gca,'YLim',[0,4.5]);
ylims = get(gca,'YLim');
plot([theta1 theta1],ylims,'--r','LineWidth',1.5);
set(gca,'YLim',ylims);
% Put a legend on it
lg1 = legend('Prior dist.','Posterior dist.','True value');
% title('Prior and posterior distributions of \theta_1');
flushLegend(lg1,'northeast');
xlabel('\theta_1');

% Second, get prior and posterior theta2
subplot(1,2,2);
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
set(gca,'YLim',[0,0.8]);
ylims = get(gca,'YLim');
plot([theta2 theta2],ylims,'--r','LineWidth',1.5);
set(gca,'YLim',ylims);
% Put a legend on it
lg2 = legend('Prior dist.','Posterior dist.','Optimal value');
% title('Prior and posterior distributions of \theta_2');
xlabel('\theta_2');
flushLegend(lg2,'northeast');
set(f,'color','white');

% Save it:
% savestr = sprintf(['FIG_discrepancy',int2str(discrep),'_results.png']);
% saveas(f,savestr);

%% Show true/optimal theta1 and theta2 locations
clc ; clearvars -except dpath ; close all ;

% Define inputs
t1min = 1.5;
t1range = 3;
t1=linspace(0,1);
t2min = 0;
t2range = 5;
% Following line gets optimal t2 for each t1 value
t2 = (4/3*(t1*t1range + t1min -1) - t2min)/t2range;

% Now get function value at optimal t2 for all t1 at x=1
fvals = dual_calib_example_fn(1,.5,.5,t1,t1min,t1range,t2,t2min,t2range,...
    0,1,0);

% Get figure showing true theta1
f = figure('pos',[20 20 650 250]);
subplot(1,2,1);
plot(t1*t1range + t1min, fvals, 'LineWidth',2);
hold on;
ylims = get(gca,'YLim');
plot([2 2],ylims,'r','LineWidth',2);
set(gca,'Ylim',ylims);
set(gca,'Xlim',[min(t1*t1range + t1min),max(t1*t1range + t1min)]);
xlabel('t_1'); ylabel('f(1,t_1,\theta_2)');
label = sprintf('True \\theta_1\nvalue');
text(2.05,0.2,label);
set(f,'color','w');

% Now get another 2D view, this time f vs theta2 at true theta1 value
subplot(1,2,2);
t2 = linspace(0,1);
fvals = dual_calib_example_fn(1,.5,.5,2,0,1,t2,t2min,t2range,...
    0,1,0);
plot(t2*t2range+t2min,fvals,'Linewidth',2);
hold on;
ylims = get(gca,'Ylim');
plot([4/3 4/3],ylims,'r','LineWidth',2);
set(gca,'Ylim',ylims);
set(gca,'Xlim',[min(t2*t2range+t2min),max(t2*t2range+t2min)]);
xlabel('t_2'); ylabel('f(1,\theta_1,t_2)');
label = sprintf('Optimal\n\\theta_2 value');
text(4/3+.1,.9,label);

% Save it
% saveas(f,'FIG_true_optimal_theta1_theta2.png');

%% Examine version without discrepancy
clc ; clearvars -except dpath ; close all ;

f=figure('pos',[20 20 700 600]);
set(f,'color','white');

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

ii=1; % for easy conversion of former code for this purpose
discrep = ii-1 ; % Select which discrepancy
Y = reshape(...
    dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,0),length(x),length(t1),length(t2));
Yd= reshape(...
    dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1,discrep),length(x),length(t1),length(t2));

% Take a look
xidx=50;
xx=reshape(X(:,xidx,:),100,100);
% Get value of x
xval = xx(1)*xrange + xmin;
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
Discrep = Yd-Y;
ea=.5;
surf(tt1*t1range+t1min,tt2*t2range+t2min,...
    reshape(Y(:,xidx,:),100,100),...
    'EdgeAlpha',ea);

%     title(sprintf('Objective function at x=%0.2f',xval));
xlabel('t_1');ylabel('t_2');zlabel('f(x,t_1,t_2)');
axis vis3d;
view([-110.4000    6.5334]);

% Save it:
% saveas(f,'FIG_obj_fn.png');

%% Results from KOH and dual calibration with various discrepancies
clc ; clearvars -except dpath ; close all ;

% Face alpha for histograms
fa = 0.6;
% Color of true/optimal value line
lcol = [218 165 32]/255 ; 

% Load results
discrep = 5; % Change this to whichever discrepancy is desired
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-21_KOH_calib_discrep',int2str(discrep)]);
load(locstr,'results');
koh_results=results;
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-15_dual_calib_discrep',int2str(discrep)]);
load(locstr,'results');


% Define inputs mins and ranges 
theta1 = 2;
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;
burn_in = results.settings.burn_in;
min_t1 = results.settings.min_t1; range_t1 = results.settings.range_t1;
min_t2 = results.settings.min_t2; range_t2 = results.settings.range_t2;

% First, get prior and posterior theta1
f = figure('pos',[10 10 600 250]);
subplot(1,2,1);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 for KOH
histogram(koh_results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',fa,'BinWidth',0.075);
% Get a histogram of theta1 for DCTO
histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',fa,'BinWidth',0.075);
% Plot true theta1
set(gca,'YLim',[0,4.5]);
ylims = get(gca,'YLim');
plot([theta1 theta1],ylims,'--','Color',lcol,'LineWidth',1.5);
set(gca,'YLim',ylims);
% Put a legend on it
lg1 = legend('Prior','KOH','DCTO','Truth');
% title('Prior and posterior distributions of \theta_1');
flushLegend(lg1,'northeast');
xlabel('t_1');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([min_t2 min_t2 + range_t2 min_t2 + range_t2 min_t2],...
    [0 0 1/range_t2 1/range_t2],'g','EdgeColor','none');
xlim([min_t2 min_t2 + range_t2]);
hold on;
% Get a histogram of theta2 for KOH
histogram(koh_results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',fa);
% Get a histogram of theta2 for DCTO
histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',fa);
% Get and plot true theta2
fmfn = @(z) dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],0,5);
set(gca,'YLim',[0,0.8]);
ylims = get(gca,'YLim');
plot([theta2 theta2],ylims,'--','Color',lcol,'LineWidth',1.5);
set(gca,'YLim',ylims);
% Put a legend on it
lg2 = legend('Prior','KOH','DCTO','Optimum');
% title('Prior and posterior distributions of \theta_2');
xlabel('t_2');
flushLegend(lg2,'northeast');
set(f,'color','white');

% Set suptitle and fix positions 
titles = {'g_0 (No discrepancy)'    ;
    'g_1 with c = 1.5'              ; 
    'g_1 with c = 3.5'              ;
    'g_2 with (c,d) = (.15,.075)'   ;
    'g_2 with (c,d) = (.65,.075)'   ;
    'g_3 with (c,d) = (.055,0)'     ;
    'g_3 with (c,d) = (.055,.1)'}   ;
% suptitle('g_3 with (c,d) = (.055,0) and informative prior on \delta_1')
suptitle(titles(discrep+1));
f.Children(2).Position(2) = .2 ; f.Children(2).Position(4) = .675 ; 
f.Children(5).Position(2) = .2 ; f.Children(5).Position(4) = .675 ; 
flushLegend(lg2,'northeast');
axes(f.Children(5)); flushLegend(lg1,'northeast');


% Save it:
savestr = ...
sprintf(['FIG_KOH_DCTO_comp_discrep',int2str(discrep),'_results.png']);
% saveas(f,savestr);
% export_fig(savestr,'-png','-m2',f);

%% Results from KOH+CTO and dual calibration with various discrepancies
clc ; clearvars -except dpath ; close all ;

% Face alpha for histograms
fa = 0.6;
% Color of true/optimal value line
lcol = [218 165 32]/255 ; 

% Load results
discrep = 4; % Change this to whichever discrepancy is desired
% Set up to load appropriate files
discrep_str = int2str(discrep);
if discrep==7 discrep = 5; discrep_str = [int2str(discrep),'_inf']; end
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_KOH_calib_discrep',discrep_str]);
load(locstr,'results');
koht1 = results.theta1 ; 
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-29_KOH_calib_step2_discrep',discrep_str]);
load(locstr,'results');
koh_results=struct('theta1',koht1,'theta2',results.theta1);
locstr = sprintf(['C:\\Users\\carle\\Documents',...
    '\\MATLAB\\NSF DEMS\\Phase 1\\',...
    'dual_calib\\dual_calib_stored_data\\'...
    '2019-03-28_DCTO_discrep',discrep_str]);
load(locstr,'results');


% Define inputs mins and ranges 
theta1 = 2;
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;
burn_in = results.settings.burn_in;
min_t1 = results.settings.min_t1; range_t1 = results.settings.range_t1;
min_t2 = results.settings.min_t2; range_t2 = results.settings.range_t2;

% First, get prior and posterior theta1
f = figure('pos',[10 10 600 250]);
subplot(1,2,1);
% Plot prior
fill([min_t1 min_t1 + range_t1 min_t1 + range_t1 min_t1],...
    [0 0 1/range_t1 1/range_t1],'g','EdgeColor','none');
xlim([min_t1 min_t1 + range_t1]);
hold on;
% Get a histogram of theta1 for KOH
histogram(koh_results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',fa,'BinWidth',0.075);
% Get a histogram of theta1 for DCTO
histogram(results.theta1(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',fa,'BinWidth',0.075);
% Plot true theta1
set(gca,'YLim',[0,4.5]);
ylims = get(gca,'YLim');
plot([theta1 theta1],ylims,'--','Color',lcol,'LineWidth',1.5);
set(gca,'YLim',ylims);
% Put a legend on it
lg1 = legend('Prior','KOH+CTO','DCTO','Truth');
% title('Prior and posterior distributions of \theta_1');
flushLegend(lg1,'northeast');
xlabel('t_1');

% Second, get prior and posterior theta2
subplot(1,2,2);
% Plot prior
fill([min_t2 min_t2 + range_t2 min_t2 + range_t2 min_t2],...
    [0 0 1/range_t2 1/range_t2],'g','EdgeColor','none');
xlim([min_t2 min_t2 + range_t2]);
hold on;
% Get a histogram of theta2 for KOH
histogram(koh_results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','r','FaceAlpha',fa);
% Get a histogram of theta2 for DCTO
histogram(results.theta2(burn_in:end,:),'Normalization','pdf',...
    'EdgeColor','none','FaceColor','b','FaceAlpha',fa);
% Get and plot true theta2
fmfn = @(z) dual_calib_example_fn(.75,0,1,theta1,0,1,z,0,1,0,1,discrep);
theta2 = fmincon(fmfn,2,[],[],[],[],0,5);
set(gca,'YLim',[0,0.9]);
ylims = get(gca,'YLim');
plot([theta2 theta2],ylims,'--','Color',lcol,'LineWidth',1.5);
set(gca,'YLim',ylims);
% Put a legend on it
lg2 = legend('Prior','KOH+CTO','DCTO','Optimum');
% title('Prior and posterior distributions of \theta_2');
xlabel('t_2');
flushLegend(lg2,'northeast');
set(f,'color','white');

% Set suptitle and fix positions 
titles = {'g_0 (No discrepancy)'    ;
    'g_1 with c = 1.5'              ; 
    'g_1 with c = 3.5'              ;
    'g_2 with (c,d) = (.15,.075)'   ;
    'g_2 with (c,d) = (.65,.075)'   ;
    'g_3 with (c,d) = (.055,0)'     ;
    'g_3 with (c,d) = (.055,.1)'}   ;
% suptitle('g_3 with (c,d) = (.055,0) and informative prior on \delta_1')
suptitle(titles(discrep+1));
f.Children(2).Position(2) = .2 ; f.Children(2).Position(4) = .675 ; 
f.Children(5).Position(2) = .2 ; f.Children(5).Position(4) = .675 ; 
flushLegend(lg2,'northeast');
axes(f.Children(5)); flushLegend(lg1,'northeast');


% Save it:
savestr = ...
sprintf(['FIG_KOHCTO_DCTO_comp_discrep',discrep_str,'_results.png']);
% saveas(f,savestr);
% export_fig(savestr,'-png','-m2',f);