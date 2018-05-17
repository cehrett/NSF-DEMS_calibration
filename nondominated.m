function nondoms = nondominated(results)
% This function takes as input a set of outcomes and returns those that are
% nondominated, ie Pareto efficient. Minimization is assumed to be optimal.

% This will hold the nondominated outcomes:
nondoms = [];

% The loop whittles down the results, keeping nondominated outcomes and
% throwing away outcomes that they dominate, until there are no results
% left.
while size(results,1) ~= 0 %

    % Get minimum on first objective
    [mi,in] = min(results(:,1));

    nondoms = [nondoms ; results(in,:) ] ; % Get new member of nondoms

    %results(in,:) = [] ; % Remove nondom from results

    % This tells us which outcomes are dominated by the lastest nondom elt.
    dominated_idx = all(results >= results(in,:),2) ;

    % Throw away the dominated outcomes
    results = results(~dominated_idx,:);

end

end