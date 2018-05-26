% Example optimization workspace
% For comparison with calibration results


% Here, I'll make a grid of costs, and then at each cost in the grid I'll
% find the curve describing the trade-off of deflection and rotation. I'll
% do this by making a grid over theta1. For fixed cost, specifying theta1
% determines theta2 as well. So I'll find deflection and rotation at each
% theta1 value in the grid. Then I'll use the results to make a plot
% describing the trade-off of deflection and rotation.
m=500;
cost_grid = linspace(15,30,m);
perfs = cell(m,1);
allperfs = [];
all_ctheta = [];

for ii = 1:m
    
    % Get theta values
    theta1 = linspace(0,(cost_grid(ii) - 15)/2,2500)'; % Restricting to 
                                                       % poss vals
    theta2 = sqrt ( (cost_grid(ii) - 15 - 2 .* theta1 ) * 4 );
    
    % for now let's just assume fixed c.
    c = 2 * ones(size(theta1));
    
    % Package inputs
    c_theta = [ c theta1 theta2 ] ;
    all_ctheta = [all_ctheta ; c_theta ] ;
    
    perf = Ex_sim(c_theta);
    
    perfs{ii}=perf;
    allperfs = [allperfs ; perf ];
    
end

% Now let's take a look at the curve for a given cost
ii=4;
plot(perfs{ii}(:,1),perfs{ii}(:,2),'-o')

% Now let's look at all costs at once
scatter3(allperfs(:,3),allperfs(:,1),allperfs(:,2));

%% Get allperfs over grid of theta1,theta2
% I don't know what I was thinking with the first grid above. I mean, I
% guess I was thinking to be able to focus on specific exact costs. But the
% calib parameter values aren't even limited to the support there! Here I
% just want a dense grid over the support: theta1 \in [0,3], theta2 \in
% [0,6]. Again let c=2. Geez!
clc; clearvars -except dpath ; close all; 
M=1000;
c=2;
theta1=linspace(0,3,M);
theta2=linspace(0,6,M);
[t1, t2] = meshgrid(theta1,theta2);
all_ctheta = [ c*ones(M^2,1) t1(:) t2(:) ] ;
allperfs = Ex_sim(all_ctheta);
ctheta_output = [all_ctheta allperfs];

% Save them here
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output'],...
    'ctheta_output');

% Now get nondominated
[nondoms, ndidx] = nondominated(allperfs);
all_ctheta_nd = all_ctheta(ndidx,:);
ctheta_output_nondom = [all_ctheta_nd nondoms ];
load([dpath,'Example\Ex_results\'...
    '2018-05-25_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');
