% Get the covariance matrix for two arrays Theta and Theta_s

function R = gp_cov (omega, X, X_s, rho, Theta, Theta_s,lambda,verbose)

% Set up matrix that will hold the covariance
R = zeros (size(X,1),size(X_s,1));

msg=0;

% If this involves only X and not Theta, then set values for rho, Theta to
% make them irrelevant:
if Theta == 0
    rho = 1; Theta = X; Theta_s = X_s ;
end


% Check whether X==X_s, Theta==Theta_s; if so, saves us half the computns
if isequal(X,X_s) && isequal(Theta,Theta_s)
    same = 1 ; % Indicator for whether X==X_2, Theta==Theta_s
else
    same = 0; 
end

for ii = 1:size(X,1)
    
    % This cuts computation in half when X,Theta == X_s,Theta_s
    if same == 1 ; termpt = ii ; else termpt = size(X_s,1) ; end

    for jj = 1:termpt
    
        R(ii,jj) = prod(omega .^(4 * ( X(ii,:) - X_s(jj,:)).^2 )) * ...
            prod (rho .^( 4 *(Theta(ii,:)-Theta_s(jj,:) ).^2 ))/lambda;
        
    end
    
    if mod(ii,100) == 0 && verbose == true
        fprintf(repmat('\b',1,msg));
        if same == 0
            msg =  fprintf('%3.1f%% done',100*ii/size(X,1));
        else
            msg = fprintf('%3.1f%% done',100*sum(1:ii)/sum(1:size(X,1)));
        end
    end
    
    
end

% If same == 1, flesh out the symmetric matrix using the lower triangle
if same == 1 ; R = (R+R') - eye(size(R,1)).*diag(diag(R)) ; end

fprintf(repmat('\b',1,msg));

end
