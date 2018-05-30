function opc = Ex_sim(c_theta,old_version)
% This function serves as an example computer model against which to assess
% the method of calibration to desired observations. Theta should be an n X
% 2 matrix, alpha and c each an n X 1 array of positive constants. The
% outputs are oscillation, performance, and cost. (We want to set
% oscillation to some constant, maximize performance, and minimize cost.)

% For backwards compatibility, there is an optional argument to use the old
% version of the function
if ~exist('old_version', 'var')
    old_version = false;
end

% Unpack input
c      = c_theta(:,1);
theta1 = c_theta(:,2);
theta2 = c_theta(:,3);

% Get outputs
% old version of oscl:
if old_version
    oscl = 1 ./ (theta1.* exp(-theta1).*abs(sin ( (theta2-1)./ c ) ) + 1) ;
else
    %new version:
    oscl = 1 ./ (theta1 .* exp(-theta1) .* ...
        exp(-abs(theta2-pi.*c/2)) + 1) ;
end

perf = 1 ./ (theta2 .^ (c - 1) .* exp ( - .75 * theta2 ) + 1) ;

cost = 15 + 2 .* theta1 + theta2 .^ 2 / 4;

% Pack up and leave
opc  = [oscl perf cost];

end