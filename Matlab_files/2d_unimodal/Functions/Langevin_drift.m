function [drift] = Langevin_drift(X,N_x,h,R,P0,m0,y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the SDE
    % N_x = dimension of the state
    % h = forward model
    % R = noise covariance
    % P0 = prior covariance
    % m0 = prior mean
    % y = observation
    %%% Output:
    % drift = drift term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute number of particles
    M = length(X)/N_x;
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % computation of sample covariances
    hx = h(Xhelp);
    xquer = mean(Xhelp,2);
    hquer = mean(hx,2);
    Pxh = 1/M*Xhelp*hx'-xquer*hquer';
    Pxx = 1/M*Xhelp*Xhelp'-xquer*xquer';
    
    % computation of the drift
    drift = Pxh*(R\(y-hx))-Pxx*(P0\(Xhelp-m0))+(N_x+1)/M*(Xhelp-xquer);
    drift = drift(:);
end

