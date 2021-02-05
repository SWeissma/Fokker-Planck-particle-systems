function [drift] = Langevin_drift_loc(X,N_x,h,R,P0,m0,y,gam)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the SDE
    % N_x = dimension of the state
    % h = forward model
    % R = noise covariance
    % P0 = prior covariance
    % m0 = prior mean
    % y = observation
    % gam = localisation parameter
    %%% Output:
    % drift = drift term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute number of particles
    M = length(X)/N_x;
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    
    hx = h(Xhelp);
    
    aprx_grad_Phi = zeros(N_x,M);
    correction_term = zeros(N_x,M);
    % computation localised drift for each particle
    for i = 1:M
        % computation of the weights
        w = exp(-1/gam*sum((Xhelp-Xhelp(:,i)).^2,1))./sum(exp(-1/gam*sum((Xhelp-Xhelp(:,i)).^2,1)));
        
        % computation of the localised sample covariances
        hx_loc = hx.*w;
        Xhelp_loc = Xhelp.*w;
        xquer_loc = sum(Xhelp_loc,2);
        gradw = w.*(Xhelp-xquer_loc);
        hxquer_loc = sum(hx_loc,2);
        Pxh_loc = w.*(Xhelp-xquer_loc)*(hx-hxquer_loc)';
        Pxx_loc = w.*(Xhelp-xquer_loc)*(Xhelp-xquer_loc)';
        
        % computation of the localised correction
        correction_term(:,i) = sum(Xhelp.*diag(Xhelp'*gradw)',2)-xquer_loc*sum(diag(Xhelp'*gradw))-sum(Xhelp.*(xquer_loc'*gradw),2)+(N_x+1)*w(i)*(Xhelp(:,i)-xquer_loc);
        % computation of the localised drift
        aprx_grad_Phi(:,i) = Pxh_loc*(R\(y-hx(:,i)))-Pxx_loc*(P0\(Xhelp(:,i)-m0));
    end

    drift = aprx_grad_Phi+correction_term;
    drift = drift(:);
end

