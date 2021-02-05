function [diff] = Langevin_diff_loc(X,N_x,gam)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the SDE
    % N_x = dimension of the state
    % gam = localisation parameter
    %%% Output:
    % diff = diffusion term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % compute number of particles
    M = length(X)/N_x;
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % computation localised diffusion for each particle
    diff = zeros(N_x*M,M);
    for i = 1:M
        w = exp(-1/gam*sum((Xhelp-Xhelp(:,i)).^2,1))./sum(exp(-1/gam*sum((Xhelp-Xhelp(:,i)).^2,1)));
        Xhelp_loc = Xhelp.*w;
        xquer_loc = sum(Xhelp_loc,2);
        Pxx_loc = w.*(Xhelp-xquer_loc)*(Xhelp-xquer_loc)';
        diff(i*N_x+1-N_x:i*N_x,i*N_x+1-N_x:i*N_x) = sqrt(2)*sqrtm(Pxx_loc);
    end
    
end

