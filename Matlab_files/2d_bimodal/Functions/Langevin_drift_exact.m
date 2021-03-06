function [drift] = Langevin_drift_exact(X,N_x,h,R,P0,m0,y,gradh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the SDE
    % N_x = dimension of the state
    % h = forward model
    % R = noise covariance
    % P0 = prior covariance
    % m0 = prior mean
    % y = observation
    % gradh = exact gradient of the forward model
    %%% Output:
    % drift = drift term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute number of particles
    M = length(X)/N_x;
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % Computation of the sample covariance
    hx = h(Xhelp);
    xquer = mean(Xhelp,2);
    Pxx = 1/M*Xhelp*Xhelp'-xquer*xquer';
    
    % Computation of the drift
    Dhx = gradh(Xhelp);
    drift = Pxx*Dhx.*(R\(y-hx))-Pxx*(P0\(Xhelp-m0))+(N_x+1)/M*(Xhelp-xquer);
    drift = drift(:);
end
