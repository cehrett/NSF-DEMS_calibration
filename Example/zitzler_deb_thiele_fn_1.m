% Zitzler–Deb–Thiele's function:
% Based on 
% https://www.researchgate.net/publication/
% 3949503_Scalable_multi-objective_optimization_test_problems
function zdt_out = zitzler_deb_thiele_fn_1(x)

    % Function assumes each observation is a row

    % Get g
    g = sum((x-.5).^2,2);
    
    % Get first output
    f1 = (1 + g) .* prod( cos( x(:,1:(end-1)) * pi/2 ),2 );
    
    % Get second output
    f2 = (1 + g) .* sin(x(:,1) * pi/2);
    
    zdt_out = [f1 f2];
    
end

% Utopia point: [0 0]
% Dystopia point: [1+.5^2*Len(x), 1+.5^2*Len(x) ]
