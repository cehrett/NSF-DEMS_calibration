% Example optimization workspace
% For comparison with calibration results


% Here, I'll make a grid of costs, and then at each cost in the grid I'll
% find the curve describing the trade-off of deflection and rotation. I'll
% do this by making a grid over theta1. For fixed cost, specifying theta1
% determines theta2 as well. So I'll find deflection and rotation at each
% theta1 value in the grid. Then I'll use the results to make a plot
% describing the trade-off of deflection and rotation.
m=25;
cost_grid = linspace(15,30,m);
perfs = cell(m,1);
allperfs = [];

for ii = 1:m
    
    % Get theta values
    theta1 = linspace(0,(cost_grid(ii) - 15)/2)'; % Restricting to poss vals
    theta2 = sqrt ( (cost_grid(ii) - 15 - 2 .* theta1 ) * 4 );
    
    % for now let's just assume fixed c.
    c = 2 * ones(size(theta1));
    
    % Package inputs
    c_theta = [ c theta1 theta2 ] ;
    
    perf = Ex_sim(c_theta);
    
    perfs{ii}=perf;
    allperfs = [allperfs ; perf ];
    
end

% Now let's take a look at the curve for a given cost
ii=4;
plot(perfs{ii}(:,1),perfs{ii}(:,2),'-o')

% Now let's look at all costs at once
scatter3(allperfs(:,3),allperfs(:,1),allperfs(:,2));
