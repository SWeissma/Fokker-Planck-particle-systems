function [ f ] = odesystem_FPPS_loc(X, h, y, R, N_x, P0, m0, B, gam)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the ODE
    % h = forward model
    % y = observation
    % R = noise covariance
    % N_x = dimension of the state
    % P0 = prior covariance
    % m0 = prior mean
    % B = Gaussian kernel covariance
    % gam = localisation parameter
    %%% Output:
    % f = drift term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute number of particles
    M = length(X)/(N_x);
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    
    % computation of the kernels and the weights
    k = zeros(M,M);
    w = zeros(M,M);
    for l = 1:M
        k(l,:) = mvnpdf(Xhelp',Xhelp(:,l)',B)';
        w(l,:) = exp(-1/gam*sum((Xhelp-Xhelp(:,l)).^2,1))./sum(exp(-1/gam*sum((Xhelp-Xhelp(:,l)).^2,1)));
    end
    
    % computation of the forward model
    hx = h(Xhelp);
    
    % computation of the localised drifts
    pre_grad_Phi = zeros(N_x,M);
    for i=1:M
        N1 = (Xhelp(:,i)-Xhelp).*k(i,:);
        hx_loc = hx.*w(i,:);
        Xhelp_loc = Xhelp.*w(i,:);
        xquer_loc = sum(Xhelp_loc,2);
        hxquer_loc = sum(hx_loc,2);
        Pxh_loc = w(i,:).*(Xhelp-xquer_loc)*(hx-hxquer_loc)';
        Pxx_loc = w(i,:).*(Xhelp-xquer_loc)*(Xhelp-xquer_loc)';
        pre_grad_Phi(:,i) = Pxh_loc*(R\(y-hx(:,i)))-Pxx_loc*(P0\(Xhelp(:,i)-m0))+Pxx_loc*((B\(mean(N1./mean(k,1),2)))+(B\(sum(N1,2)/sum(k(i,:)))));
    end
    
    f = pre_grad_Phi(:);        
    
end