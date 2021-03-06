function [ f ] = odesystem_FPPS(X, h, y, R, N_x, P0, m0, B)
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
    % gradh = exact gradient of the forward model
    %%% Output:
    % f = drift term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M = length(X)/(N_x);
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % computation of the kernels 
    k = zeros(M,M);
    for j = 1:M
        k(j,:) = mvnpdf(Xhelp',Xhelp(:,j)',B)';
    end
    
    % computation of the sample covariances
    hx = h(Xhelp);
    hxquer = mean(hx,2);
    uquer = mean(Xhelp,2);
    Pxx = 1/M*Xhelp*Xhelp'-uquer*uquer';
    Pxh = 1/M*Xhelp*hx'-uquer*hxquer';
    
    
    % computation of the drift
    K1 = zeros(N_x,M);
    K2 = zeros(N_x,M);
    for i=1:M
        N1 = (Xhelp(:,i)-Xhelp).*k(i,:);
        K1(:,i) = sum(N1,2)/sum(k(i,:));
        K2(:,i) = mean(N1./mean(k,1),2);
    end  
    f = (Pxh*(R\(y-hx))-Pxx*(P0\(Xhelp-m0)))+Pxx*(B\K1)+Pxx*(B\K2);    
    f = f(:);
    
    
end