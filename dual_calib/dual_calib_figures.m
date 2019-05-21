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

% f=figure('pos',[640 5 540 800]);
% set(f,'color','white');

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
%     subplot(3,2,ii);
    if mod(ii,2)==1 figure('pos',[10 + ii*20, 5, 540, 250],'color','w');end
    subplot(1,2,mod(ii-1,2)+1);
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
    zlim([0,1.33]);
    
    % Sort out title and labels
    dtc = discrep_title_content(:,ii);
    title(sprintf('g_%d, %s',dtc{:}));
    xlabel('t_1');ylabel('t_2');zlabel('f(x,t_1,t_2)');
    axis vis3d;
    view([-110.4000    6.5334]);
    
    % Save
    savestr = sprintf(['FIG_obj_fn_g',int2str(ceil(ii/2))]);
    if mod(ii,2)==0 export_fig(savestr,'-png','-m2'); end
end

%%% Fix sizing
% f.Children(6).Position = [0.0 0.685 0.575 0.32];
% f.Children(5).Position = [0.5 0.685 0.575 0.32];
% f.Children(4).Position = [0.0 0.355 0.575 0.32];
% f.Children(3).Position = [0.5 0.355 0.575 0.32];
% f.Children(2).Position = [0.0 0.025 0.575 0.32];
% f.Children(1).Position = [0.5 0.025 0.575 0.32];

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

%% Get histograms of posterior predictive distributions
clc ; clearvars -except dpath ; close all ;

% Set real theta1, optimal theta2, discrepancy version, whether modular,
% and known observation variance
sig2=0.05;
true_theta1 = 2;
% discrep = 1;
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

% Set desired observations
n = 3 ; % Number of points to use for integration
x = linspace(0,1,n)'; % Get points 


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

% Add discrepancy means to computer model output, rescale sds
discrep_gp_post_sds = discrep_gp_post_sds_std * std_y ;% Put on orig scale
discrep_gp_post_means = discrep_gp_post_means_std * std_y;% Put on orig scl
posterior_preds = comp_model_output + discrep_gp_post_means;
posterior_lb = posterior_preds - 2*discrep_gp_post_sds; % lower bound
posterior_ub = posterior_preds + 2*discrep_gp_post_sds; % upper bound

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
f=figure('pos',[10 10 800 250]);
for ii = 1:n
    posterior_distro = ...
        mean(normpdf(linspace(0,1),...
        posterior_preds(:,ii),discrep_gp_post_sds(:,ii)));
    subplot(1,n,ii);
%     histogram(prior_model_output(:,ii),'Normalization','pdf',...
%         'DisplayStyle','stairs','LineWidth',2); hold on;
    histogram(prior_model_output(:,ii),'Normalization','pdf',...
        'EdgeColor','none','FaceColor','b','FaceAlpha',.5); hold on;
    histogram(prior_model_output(:,ii),'Normalization','pdf',...
        'DisplayStyle','stairs','EdgeColor','k');
%     plot(linspace(0,1),posterior_distro,'LineWidth',2); 
    area(linspace(0,1),posterior_distro,'FaceColor','r','FaceAlpha',.5); 
    ylims = get(gca,'YLim');
    plot(true_output(ii)*[1 1],ylims,'--k','LineWidth',1.5);
    title(sprintf('x = %g',x(ii)*xrange+xmin));
    xlim([0,1]);
    set(gca,'YTick',[]);
    set(gca,'YLim',ylims);
end
% Set title
if discrep_str == '5_inf' discrep=7 ; end % To get proper title
titles = {'g_0 (No discrepancy)'    ;
    'g_1 with c = 1.5'              ;
    'g_1 with c = 3.5'              ;
    'g_2 with (c,d) = (.15,.075)'   ;
    'g_2 with (c,d) = (.65,.075)'   ;
    'g_3 with (c,d) = (.055,0)'     ;
    'g_3 with (c,d) = (.055,.1)'    ;
    'g_3 with (c,d) = (.055,0) and inf. prior'}     ;
suptitle(titles(discrep+1));


% Save it:
set(f,'Color','w');
savestr = ...
sprintf(['FIG_DCTO_post_pred_dist_discrep',discrep_str,'.png']);
% saveas(f,savestr);
export_fig(savestr,'-png','-m2',f);