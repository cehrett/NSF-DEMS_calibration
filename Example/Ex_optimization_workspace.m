% Example optimization workspace
% For comparison with calibration results


% Here, I'll make a grid of costs, and then at each cost in the grid I'll
% find the curve describing the trade-off of deflection and rotation. I'll
% do this by making a grid over theta1. For fixed cost, specifying theta1
% determines theta2 as well. So I'll find deflection and rotation at each
% theta1 value in the grid. Then I'll use the results to make a plot
% describing the trade-off of deflection and rotation.
m=8;
cost_grid = linspace(15,30,m);

for ii = 1:m
    
    % Get theta values
    theta1 = linspace(0,(cost_grid(ii) - 15)/2);
    theta2 = sqrt ( (cost_grid(ii) - 15 - 2 .* theta1 ) * 4 );
    
end