function [diff] = Langevin_diff(X,N_x)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the SDE
    % N_x = dimension of the state
    %%% Output:
    % diff = diffusion term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % compute number of particles
    M = length(X)/N_x;
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % computation of the sample covariance
    xquer = mean(Xhelp,2);
    Pxx = 1/M*Xhelp*Xhelp'-xquer*xquer';
    
    % computation of the diffusion
    diff_h = sqrt(2)*sqrtm(Pxx);
    diff = kron(eye(M),diff_h);
end
