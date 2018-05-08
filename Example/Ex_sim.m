function [oscl perf cost] = Ex_sim(theta, alpha, c)
% This function serves as an example computer model against which to assess
% the method of calibration to desired observations. Theta should be an n X
% 2 matrix, alpha and c each an n X 1 array of positive constants. The
% outputs are oscillation, performance, and cost. (We want to set
% oscillation to some constant, maximize performance, and minimize cost.)

theta1=theta(:,1);
theta2=theta(:,2);

oscl = theta1 .* sin ( theta2 ./ c ) ;

perf = theta2 .^ (alpha - 1) .* exp ( - theta2 ) ;

cost = theta1 + theta2 .^ 2 ;

end