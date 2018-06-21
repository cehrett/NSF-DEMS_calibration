% Get the covariance matrix for two arrays Theta and Theta_s

function R = gp_cov(omega, X, X_s, rho, Theta, Theta_s,lambda,verbose)

% If this involves only X and not Theta, then set values for rho, Theta to
% make them irrelevant:
if Theta == 0
    rho = 1; Theta = X; Theta_s = X_s ;
end

% Initialize the covariance matrix
R = ones(size(X,1),size(X_s,1));

% % Get matrices of X, X_s for vectorized computations
% BX = repmat(X,1,size(X_s,1));
% BX_s = repmat(X_s',size(X,1),1);
% ExpX = 4 * (BX-BX_s).^2;
% % Now same for Theta
% BTheta = repmat(Theta,1,size(Theta_s,1));
% BTheta_s = repmat(Theta_s',size(Theta,1),1);
% ExpT = 4 * (BTheta - BTheta_s).^2;

msg=0; % used for verbose output
% Now loop through all the elements of omega.
for ii = 1:numel(omega)
    BX = repmat(X(:,ii),1,size(X_s,1));
    BX_s = repmat(X_s(:,ii)',size(X,1),1);
    ExpX = 4 * (BX-BX_s).^2;
    Omega = omega(ii) * ones(size(R)); % Get a matrix of the current omega
    R = R .* Omega.^ExpX ; 
    if verbose
        percdone = 100*ii/sum([numel(omega),numel(rho)]);
        fprintf(repmat('\b',1,msg));
        msg=fprintf('%3.2f%% done\n',percdone);
    end
end
% Now do the same for all elements of rho.
for ii = 1:numel(rho)
    BTheta = repmat(Theta(:,ii),1,size(Theta_s,1));
    BTheta_s = repmat(Theta_s(:,ii)',size(Theta,1),1);
    ExpT = 4 * (BTheta - BTheta_s).^2;
    Rho = rho(ii) * ones(size(R));
    R = R .* Rho.^ExpT ;
    if verbose
        percdone = 100*(ii+numel(omega))/sum([numel(omega),numel(rho)]);
        fprintf(repmat('\b',1,msg));
        msg=fprintf('%3.2f%% done\n',percdone);
    end
end

% Use marginal precision:
R = R/lambda;

end
