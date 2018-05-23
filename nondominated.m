function [nondoms,ndidx] = nondominated(results)
% This function takes as input a set of outcomes and returns those that are
% nondominated, ie Pareto efficient. Minimization is assumed to be optimal.

% This will hold the nondominated outcomes:
nondoms = [];

% This will hold the indices of those outcomes:
ndidx = [];

% The loop whittles down the results, keeping nondominated outcomes and
% throwing away outcomes that they dominate, until there are no results
% left.
still_left = 1 ; % Will be set to zero to end loop
while still_left > 0 %

    % Get minimum on first objective
    [mi,in] = min(results(:,1));

    nondoms = [nondoms ; results(in,:) ] ; % Get new member of nondoms
    ndidx   = [ndidx   ; in            ] ;

    % This tells us which outcomes are dominated by the lastest nondom elt.
    dominated_idx = all(results >= results(in,:),2) ;

    % Throw away the dominated outcomes
    results(dominated_idx,:)=Inf*ones(sum(dominated_idx),size(results,2));
    
    still_left = sum(results(:,1) ~= Inf);
    
    if rand < 0.01
        fprintf('Still left: %d\n',still_left);
    end

end

end