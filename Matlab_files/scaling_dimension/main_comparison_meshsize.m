%%% Fokker--Planck Particle System: Toy example for high dimension

clear;
close all;

% path to functions
path('./Functions',path);


% options for the ode setup
tspan1 = [0 1000]; % run ODEs until time T=1000
tspan2 = [0 1000];

options = odeset('RelTol',1e-8,'AbsTol',1e-8);  % accuracy of ODE solver


% Properties of the figures
fig1 = figure(1);
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig2 = figure(2);
clf(fig2)
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig3 = figure(3);
clf(fig3)
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig4 = figure(4);
clf(fig4)
set(fig4, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig5 = figure(5);
clf(fig5)
set(fig5, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig6 = figure(6);
clf(fig6)
set(fig6, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig7 = figure(7);
clf(fig7)
set(fig7, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig8 = figure(8);
clf(fig8)
set(fig8, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);
fig9 = figure(9);
clf(fig9)
set(fig9, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.8]);


Dim = [4,6,8];  % considered Dimensions (discretization of the PDE)
MM = [50,100,200];  % considered ensemble sizes
N_x = 4;    % truncation of the KL expansion

mean_m1 = zeros(length(Dim),length(MM));     % estimated mean for B1
trace_P1 = zeros(length(Dim),length(MM));    % estimated trace of covariance for B1 

mean_m2 = zeros(length(Dim),length(MM));    % estimated mean for B2
trace_P2 = zeros(length(Dim),length(MM));   % estimated trace of covariance for B2

mean_m3 = zeros(length(Dim),length(MM));    % estimated mean for B3
trace_P3 = zeros(length(Dim),length(MM));   % estimated trace of covariance for B3

trace_theoretical = zeros(length(Dim),1);
i = 1; % running variable for subplots
n = 1; % running variable for subplots


D = zeros(N_x,N_x);
for j = 1:N_x
    D(j,j) = (j^2)^(-1);
end
m_ref = sqrt(D)*randn(N_x,1);

m0 = zeros(N_x,1);                % prior mean
P0 = 1*eye(N_x,N_x)*D;     % prior covariance


for ell = Dim % dimension of target distribution
    m=1;
    for M = MM % ensemble size
        
        l = ell;          % discretization of the PDE
        
        I = 2^l;            % dimension of the discretized system
        xx_b = linspace(0,1,I+1);       % domain with boundary
        xx = xx_b(2:end-1);             % domain without boundary
        
        
        % KL Basis
        V = zeros(I-1,N_x);
        for j = 1:N_x
            V(:,j) = sqrt(2*pi)*sin(j*pi*xx)';
        end
        
        % evaluation of the KL expansion
        u = @(xi) V*xi;
        
        % evaluate reference solution
        y = u(m_ref);
        R = 10*eye(I-1,I-1);
        
        % theoretical posterior mean and covariance
        m_ast = m0 + P0*V'*((V*P0*V'+R)\(y-V*m0));
        P_ast = P0 - P0*V'*((V*P0*V'+R)\V)*P0;

        % sample from posterior
        sample_size = 10000;
        sample_posterior = m_ast + chol(P_ast)'*randn(N_x,sample_size);
        
        trace_theoretical(n,1) = trace(P_ast);
        
        u_posterior_est = u(sample_posterior);
        u_posterior_est_mean = mean(u_posterior_est,2);

        SD_u_posterior = sqrt(var(u_posterior_est')');
        LB_u_posterior = mean(u_posterior_est,2)-SD_u_posterior;
        UB_u_posterior = mean(u_posterior_est,2)+SD_u_posterior;

        
        
        % prior sample
        x0 = mvnrnd(m0,P0,M)';  % initial ensemble

        % optimal kernel bandwidth wrt. AMISE
        delta = 1/(N_x+4);
        kernel_width = (4/(N_x+2))^delta*(1/M)^delta; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% Kernel choice 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % scaling of the Gaussian kernel
        B1 = kernel_width*P0;
        
        %%%%%%%%%%%%%%%%%%%% run the algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [t1, X1] = ode45(@(t,X) odesystem_FPPS_fix_u(X,u,N_x,R,y,m0,P0,B1), tspan1, x0(:),options);
        x1 = reshape(X1(end,:),[N_x,M]);
        
        runningtime1 = toc
        
        %%%%%%%%%%%%%%%%%%%% evaluate the algorithm %%%%%%%%%%%%%%%%%%%%%%%
        
        X_est1 = x1;
        u_est1 = u(X_est1);
        x_quer = mean(x1,2);
        P_t = 1/M*(x1-x_quer)*(x1-x_quer)'; % sample covariance
        
        % Resampling to produce approximate sample of the posterior
        Z = sampling_kernel(X_est1,B1,1000);   % sampling from the kernels
        X_sample1 = resampling(Z,X_est1,B1,10000);    % resampling
        
        % evaluate KL expansion
        u_1_est = u(X_sample1);
        u_1_est_mean = mean(u_1_est,2);
        
        % compute empirical SD
        SD_u1 = sqrt(var(u_1_est')');
        LB_u1 = mean(u_1_est,2)-SD_u1;
        UB_u1 = mean(u_1_est,2)+SD_u1;
        
        est_P = 1/10000*X_sample1*X_sample1'-mean(X_sample1,2)*mean(X_sample1,2)';    % covariance estimate
        
        % trace of the sample covariance
        trace_P1(n,m) = trace(est_P);
        
        % estimated mean
        mean_m1(n,m) = mean(mean(X_sample1,2));   
        
        % plot only a sample size of 512
        X_sample1 = X_sample1(:,1:512);
        
        %%%%%%%%%%%%%%%%%%%% numerical result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % compute the kernel
        k = zeros(M,M);
        for j = 1:M
            k(j,:) = 1/sqrt(det(B1))*exp(-1/2*sum((sqrtm(B1)\(X_est1-X_est1(:,j))).^2,1));
        end
    
        K1 = zeros(N_x,M);
        K2 = zeros(N_x,M);
        for l=1:M
            N1 = (B1\(X_est1(:,l)-X_est1)).*k(l,:);
            K1(:,l) = sum(N1,2)/sum(k(l,:));
            K2(:,l) = sum(N1./sum(k,1),2);
        end
        uquer = mean(u(X_est1),2);
        Pu_t = 1/M*(X_est1-x_quer)*(u(X_est1)-uquer)';
        
        % force which pushes the particle system into MAP direction
        grad_V_1 = -(Pu_t*(R\(u(X_est1)-y)))-P_t*(P0\(X_est1-m0));
        % force which spreads the particles
        grad_V_2 = (P_t*(K1)+P_t*(K2));
        
        %%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(1)
        subplot(length(Dim),length(MM),i)
        scatter(X_est1(1,:),-grad_V_1(1,:),'MarkerEdgeColor',[0 0 0.3],'DisplayName','gradient');hold on
        scatter(X_est1(1,:),grad_V_2(1,:),'x','MarkerFaceColor',[0 1 0],'DisplayName','gradient','LineWidth',1);hold on

        str = sprintf('$B, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
         
        figure(2)
        subplot(length(Dim),length(MM),i)
        scatter(sample_posterior(1,:),sample_posterior(2,:),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','posterior');hold on
        s2 = scatter(X_sample1(1,:),X_sample1(2,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
        scatter(X_est1(1,:),X_est1(2,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$B=B_1$');hold on
        ylabel('x_{2}')
        xlabel('x_1')
        str = sprintf('$B, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        l = legend('show','Location','southwest');
        set(l,'Interpreter','latex','FontSize',14);


        figure(7)
        subplot(length(Dim),length(MM),i)
        Uest_sample_mean = mean(u_posterior_est,2);

        SD_sample = sqrt(var(u_posterior_est')');
        LB_sample = mean(u_posterior_est,2)-SD_sample;
        UB_sample = mean(u_posterior_est,2)+SD_sample;
        p1 = plot(xx_b,[0;Uest_sample_mean;0],'-','Color',[0 0 0.5],'LineWidth',3,'DisplayName','estimate');hold on
        plot(xx_b,[0;LB_sample;0],'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
        plot(xx_b,[0;UB_sample;0],'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
        fill([xx;xx],[LB_sample';UB_sample'],'--','edgecolor',[0 0 1],'LineWidth',0.5);
        fill([xx',xx']',[LB_sample,UB_sample]','--','edgecolor',[0 0 1],'LineWidth',0.5);
        
        
        p2 = plot(xx_b,[0;u_1_est_mean;0],'-.','Color',[0.5 0 0],'LineWidth',3,'DisplayName','posterior');hold on
        plot(xx_b,[0;LB_u1;0],':','color',[1 0 0],'LineWidth',3,'DisplayName','credible set');hold on
        plot(xx_b,[0;UB_u1;0],':','color',[1 0 0],'LineWidth',3,'DisplayName','credible set');hold on
        fill([xx;xx],[LB_u1';UB_u1'],':','edgecolor',[1 0 0],'LineWidth',0.5);
        fill([xx',xx']',[LB_u1,UB_u1]',':','edgecolor',[1 0 0],'LineWidth',0.5);

        
        str = sprintf('$B, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        xlabel('x')
        legend('show',[p1,p2],'FontSize',16,'Location','southeast')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% Kernel choice 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % kernel scaling
        B2 = kernel_width*diag(diag(P_ast));
        
        %%%%%%%%%%%%%%%%%%%% run the algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [t2, X2] = ode45(@(t,X) odesystem_FPPS_fix_u(X,u,N_x,R,y,m0,P0,B2), tspan1, x0(:),options);
        x2 = reshape(X2(end,:),[N_x,M]);
        runningtime2 = toc
        
        %%%%%%%%%%%%%%%%%%%% evaluate the algorithm %%%%%%%%%%%%%%%%%%%%%%%
        X_est2 = x2;
        u_est2 = u(X_est2);
        x_quer = mean(x2,2);
        P_t = 1/M*(x2-x_quer)*(x2-x_quer)'; % sample covariance

        % Resampling to produce approximate sample of the posterior
        Z = sampling_kernel(X_est2,B2,1000);    % sampling from the kernels
        X_sample2 = resampling(Z,X_est2,B2,10000);    % resampling
        
        % evaluate KL expansion
        u_2_est = u(X_sample2);
        u_2_est_mean = mean(u_2_est,2);
        
        % compute empirical SD
        SD_u2 = sqrt(var(u_2_est')');
        LB_u2 = mean(u_2_est,2)-SD_u2;
        UB_u2 = mean(u_2_est,2)+SD_u2;
        
        est_P = 1/10000*X_sample2*X_sample2'-mean(X_sample2,2)*mean(X_sample2,2)';
        
        % trace of the sample covariance
        trace_P2(n,m) = trace(est_P);
        
        % estimated mean
        mean_m2(n,m) = mean(mean(u_2_est,2));
        
        X_sample2 = X_sample2(:,1:512);

        %%%%%%%%%%%%%%%%%%%% numerical result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % compute the kernel
        k = zeros(M,M);
        for j = 1:M
            k(j,:) = 1/sqrt(det(B2))*exp(-1/2*sum((sqrtm(B2)\(X_est2-X_est2(:,j))).^2,1));
        end
    
        K1 = zeros(N_x,M);
        K2 = zeros(N_x,M);
        for l=1:M
            N1 = (B2\(X_est2(:,l)-X_est2)).*k(l,:);
            K1(:,l) = sum(N1,2)/sum(k(l,:));
            K2(:,l) = sum(N1./sum(k,1),2);
        end
        uquer = mean(u(X_est2),2);
        Pu_t = 1/M*(X_est2-x_quer)*(u(X_est2)-uquer)';
        % force which pushes the particle system into MAP direction
        grad_V_1 = -(Pu_t*(R\(u(X_est2)-y)))-P_t*(P0\(X_est2-m0));
        % force which spreads the particles
        grad_V_2 = (P_t*(K1)+P_t*(K2));

        
        %%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(3)
        subplot(length(Dim),length(MM),i)
        scatter(X_est2(1,:),-grad_V_1(1,:),'MarkerEdgeColor',[0 0 0.3],'DisplayName','gradient');hold on
        scatter(X_est2(1,:),grad_V_2(1,:),'x','MarkerFaceColor',[0 1 0],'DisplayName','gradient','LineWidth',1);hold on

        str = sprintf('$B_*, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)

        
        figure(4)
        subplot(length(Dim),length(MM),i)
        scatter(sample_posterior(1,:),sample_posterior(2,:),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','posterior');hold on
        scatter(X_sample2(1,:),X_sample2(2,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
        scatter(X_est2(1,:),X_est2(2,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$B=B_2$');hold on
        ylabel('x_{2}')
        xlabel('x_1')
        str = sprintf('$B_*, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        l = legend('show','Location','southwest');
        set(l,'Interpreter','latex','FontSize',14);
        
        
        figure(8)
        subplot(length(Dim),length(MM),i)
        Uest_sample_mean = mean(u_posterior_est,2);

        SD_sample = sqrt(var(u_posterior_est')');
        LB_sample = mean(u_posterior_est,2)-SD_sample;
        UB_sample = mean(u_posterior_est,2)+SD_sample;
        p1 = plot(xx_b,[0;Uest_sample_mean;0],'-','Color',[0 0 0.5],'LineWidth',3,'DisplayName','posterior');hold on
        plot(xx_b,[0;LB_sample;0],'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
        plot(xx_b,[0;UB_sample;0],'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
        fill([xx;xx],[LB_sample';UB_sample'],'--','edgecolor',[0 0 1],'LineWidth',0.5);
        fill([xx',xx']',[LB_sample,UB_sample]','--','edgecolor',[0 0 1],'LineWidth',1);
        
        
        p2 = plot(xx_b,[0;u_2_est_mean;0],'-.','Color',[0.5 0 0],'LineWidth',3,'DisplayName','estimate');hold on
        plot(xx_b,[0;LB_u2;0],':','color',[1 0 0],'LineWidth',3,'DisplayName','credible set');hold on
        plot(xx_b,[0;UB_u2;0],':','color',[1 0 0],'LineWidth',3,'DisplayName','credible set');hold on
        fill([xx;xx],[LB_u2';UB_u2'],':','edgecolor',[1 0 0],'LineWidth',0.5);
        fill([xx',xx']',[LB_u2,UB_u2]',':','edgecolor',[1 0 0],'LineWidth',0.5);
        
        str = sprintf('$B_*, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        xlabel('x')
        legend('show',[p1,p2],'FontSize',16,'Location','southeast')
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% Kernel choice 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % kernel scaling applied directly in the ode function
        
        %%%%%%%%%%%%%%%%%%%% run the algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [t3, X3] = ode45(@(t,X) odesystem_FPPS_u(X,u,N_x,R,y,m0,P0,delta), tspan2, x0(:),options);
        x3 = reshape(X3(end,:),[N_x,M]);
        runningtime3 = toc

        %%%%%%%%%%%%%%%%%%%% evaluate the algorithm %%%%%%%%%%%%%%%%%%%%%%%
        X_est3 = x3;
        u_est3 = u(X_est3);
        x_quer = mean(x3,2);
        P_t = 1/M*(x3-x_quer)*(x3-x_quer)'; % sample covariance
        B_est = diag(P_t);
        B3 = kernel_width*(diag(B_est));
        
        
        % Resampling to produce approximate sample of the posterior
        Z = sampling_kernel(X_est3,B3,1000);   % sampling from the kernels
        X_sample3 = resampling(Z,X_est3,B3,10000);    % resampling
        
        % evaluate KL expansion
        u_3_est =u(X_sample3);
        u_3_est_mean = mean(u_3_est,2);

        % compute empirical SD
        SD_u3 = sqrt(var(u_3_est')');
        LB_u3 = mean(u_3_est,2)-SD_u3;
        UB_u3 = mean(u_3_est,2)+SD_u3;
        
        est_P = 1/10000*X_sample3*X_sample3'-mean(X_sample3,2)*mean(X_sample3,2)';    % covariance estimate
        
        % trace of the sample covariance
        trace_P3(n,m) = trace(est_P);
        
        % estimated mean
        mean_m3(n,m) = mean(mean(u_3_est,2));
        
        % plot only a sample size of 512
        X_sample3 = X_sample3(:,1:512);
        
        %%%%%%%%%%%%%%%%%%%% numerical result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute the kernel
        k = zeros(M,M);
        for j = 1:M
            k(j,:) = 1/sqrt(det(B3))*exp(-1/2*sum((sqrtm(B3)\(X_est3-X_est3(:,j))).^2,1));
        end

        K1 = zeros(N_x,M);
        K2 = zeros(N_x,M);
        for l=1:M
            N1 = (B3\(X_est3(:,l)-X_est3)).*k(l,:);
            K1(:,l) = sum(N1,2)/sum(k(l,:));
            K2(:,l) = sum(N1./sum(k,1),2);
        end
        
        uquer = mean(u(X_est3),2);
        Pu_t = 1/M*(X_est3-x_quer)*(u(X_est3)-uquer)';
        % force which pushes the particle system into MAP direction
        grad_V_1 = -(Pu_t*(R\(u(X_est3)-y)))-P_t*(P0\(X_est3-m0));
        % force which spreads the particles
        grad_V_2 = (P_t*(K1)+P_t*(K2));
        
        %%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(5)
        subplot(length(Dim),length(MM),i)
        scatter(X_est3(1,:),-grad_V_1(1,:),'MarkerEdgeColor',[0 0 0.3],'DisplayName','gradient');hold on
        scatter(X_est3(1,:),grad_V_2(1,:),'x','MarkerFaceColor',[0 1 0],'DisplayName','gradient','LineWidth',1);hold on

        str = sprintf('$B_t, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        
        figure(6)
        subplot(length(Dim),length(MM),i)
        scatter(sample_posterior(1,:),sample_posterior(2,:),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','posterior');hold on
        scatter(X_sample3(1,:),X_sample3(2,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
        scatter(X_est3(1,:),X_est3(2,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$B=B_3$');hold on
        
        ylabel('x_{2}')
        xlabel('x_1')
        str = sprintf('$B_t, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        l = legend('show','Location','southwest');
        set(l,'Interpreter','latex','FontSize',14);

        figure(9)
        subplot(length(Dim),length(MM),i)
        
        Uest_sample_mean = mean(u_posterior_est,2);
        
        SD_sample = sqrt(var(u_posterior_est')');
        LB_sample = mean(u_posterior_est,2)-SD_sample;
        UB_sample = mean(u_posterior_est,2)+SD_sample;
        p1 = plot(xx_b,[0;Uest_sample_mean;0],'-','Color',[0 0 0.5],'LineWidth',3,'DisplayName','estimate');hold on
        plot(xx_b,[0;LB_sample;0],'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
        plot(xx_b,[0;UB_sample;0],'--','color',[0 0 1],'LineWidth',3,'DisplayName','credible set');hold on
        fill([xx;xx],[LB_sample';UB_sample'],'--','edgecolor',[0 0 1],'LineWidth',0.5);
        fill([xx',xx']',[LB_sample,UB_sample]','--','edgecolor',[0 0 1],'LineWidth',1);
        
        
        p2 = plot(xx_b,[0;u_3_est_mean;0],'-.','Color',[0.5 0 0],'LineWidth',3,'DisplayName','posterior');hold on
        plot(xx_b,[0;LB_u3;0],':','color',[1 0 0],'LineWidth',3,'DisplayName','credible set');hold on
        plot(xx_b,[0;UB_u3;0],':','color',[1 0 0],'LineWidth',3,'DisplayName','credible set');hold on
        fill([xx;xx],[LB_u3';UB_u3'],':','edgecolor',[1 0 0],'LineWidth',0.5);
        fill([xx',xx']',[LB_u3,UB_u3]',':','edgecolor',[1 0 0],'LineWidth',0.5);
        
        str = sprintf('$B_t, N_y=%d, M=%d$',2^ell,M);
        title(str,'Interpreter','latex','FontSize',20)
        xlabel('x')
        legend('show',[p1,p2],'FontSize',16,'Location','southeast')

        m = m+1;
        i= i+1;
    end
    n = n+1;
    
end
  
