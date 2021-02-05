function [X_path,time] = Langevin_SMC(h,P0,x0,y,T,R,scale)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % h = forward model
    % P0 = prior covariance
    % x0 = initial ensemble
    % y = observation
    % T = final time
    % R = noise covariance
    % scale = scaling of the loglikelihood for SMC
    %%% Output:
    % X_path = (approximated) path of the SDE
    % time = corresponding time nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Xhelp = x0;
    [N_x,M] = size(x0);

    time = zeros(1,1);
    X_path = zeros(1,N_x*M);
    
    PP = eye(M)-ones(M)/M;
    P0_inv = P0\eye(N_x,N_x);
    K = length(y);
    R_inv = R\eye(K,K);
    
    current_time = 0;
    j=1;
    while current_time < T % run Langevin until time T
        hx = h(Xhelp);
        xquer = mean(Xhelp,2);
        A = 1/sqrt(M)*Xhelp*PP;
        B = 1/sqrt(M)*hx*PP;
        
        % compute the drift
        drift = -scale*A*B'*(R_inv*(hx-y))-A*A'*P0_inv*Xhelp+(N_x+1)/M*(Xhelp-xquer);
        
        % compute adaptive step size 
        max_drift = max(sqrt(sum(drift.^2,2)));
        delta_t = 1/max_drift;
        if (current_time + delta_t) >= T
            delta_t = T-current_time;
        end
        
        % compute Euler-Maruyama update step
        Xhelp = Xhelp + delta_t*drift+ sqrt(2*delta_t)*A*randn(M,M);
        
        % save the path
        X_path(j,:) = Xhelp(:);
        
        % update current time
        current_time = current_time+delta_t;
        
        % save time
        time(j) = current_time;
        
        
        j = j+1;
    end    
    
end



