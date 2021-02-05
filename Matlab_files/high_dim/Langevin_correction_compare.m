function [X_path,time] = Langevin_correction_compare(h,P0,x0,y,T,R,TT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % h = forward model
    % P0 = prior covariance
    % x0 = initial ensemble
    % y = observation
    % T = final time
    % R = noise covariance
    % TT = time points to hit for simulation
    %%% Output:
    % X_path = (approximated) path of the SDE
    % time = corresponding time nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize the state
    Xhelp = x0;
    
    % get the dimension and ensemble size
    [N_x,M] = size(x0);

    % initialize time and path for the SDE
    time = zeros(1,1);
    X_path = zeros(1,N_x*M);
    
    % auxilary object for computation of the sample covariance
    PP = eye(M)-ones(M)/M;
    
    % compute inverse of prior and noise covariances
    P0_inv = P0\eye(N_x,N_x);
    K = length(y);
    R_inv = R\eye(K,K);
    
    % initialize current time
    current_time = 0;
    j=1;
    count = 1;
    
    % run the Langevin dynamics up to final time T
    while current_time < T
        
        % evaluate the forward model
        hu = h(Xhelp);
        % ensemble mean
        uquer = mean(Xhelp,2);
        
        % auxilary objects for computation of sample covariance
        A = 1/sqrt(M)*Xhelp*PP;
        B = 1/sqrt(M)*hu*PP;
        
        % compute the drift
        drift = -A*B'*(R_inv*(hu-y))-A*A'*P0_inv*Xhelp+(N_x+1)/M*(Xhelp-uquer);
        
        % compute adaptive step size
        max_drift = max(sqrt(sum(drift.^2,2)));
        delta_t = 1/max_drift;
        
        % hit the times in TT
        if (current_time + delta_t) >= TT(count)
            delta_t = TT(count)-current_time;
            count = count+1;
        end
        
        % update by Euler-Maruyama method
        Xhelp = Xhelp + delta_t*drift+ sqrt(2*delta_t)*A*randn(M,M);%sqrtm(Cuu)*randn(I,J);
        
        % save the current iteration as path
        X_path(j,:) = Xhelp(:);
        
        % update time
        current_time = current_time+delta_t;
        
        % save the time
        time(j) = current_time;
        
        j = j+1;
    end    
    
end
