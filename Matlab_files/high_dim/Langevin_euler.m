function [X_path,time] = Langevin_euler(h,P0,x0,y,T,R)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % h = forward model
    % P0 = prior covariance
    % x0 = initial ensemble
    % y = observation
    % T = final time
    % R = noise covariance
    %%% Output:
    % X_path = (approximated) path of the SDE
    % time = corresponding time nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Xhelp = x0;
    [I,J] = size(x0);

    time = zeros(1,1);
    X_path = zeros(1,I*J);
    
    PP = eye(J)-ones(J)/J;
    P0_inv = P0\eye(I,I);
    K = length(y);
    R_inv = R\eye(K,K);
    
    current_time = 0;
    j=1;
    while current_time < T
        hx = h(Xhelp);

        A = 1/sqrt(J)*Xhelp*PP;
        B = 1/sqrt(J)*hx*PP;

        drift = -A*B'*(R_inv*(hx-y))-A*A'*P0_inv*Xhelp;
        max_drift = max(sqrt(sum(drift.^2,2)));
        delta_t = 0.05/max_drift;
        if (current_time + delta_t) >= T
            delta_t = T-current_time;
        end
        
        Xhelp = Xhelp + delta_t*drift+ sqrt(2*delta_t)*A*randn(J,J);
        X_path(j,:) = Xhelp(:);
        current_time = current_time+delta_t;
        time(j) = current_time;
        
        
        j = j+1;
    end    
end

