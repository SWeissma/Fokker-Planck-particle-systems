%%% Evaluation FPPS
clear;
close all;

% load the corresponding data:
% note that one has to run FPPS.m and for M2 = 128,256,512 first.
load('reference_data.mat');
load('MCMC.mat')

%%% scatter plots

fig2 = figure(2);
clf(fig2)
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.9]);

load('FPPS_J128.mat');
Xi_est = reshape(Xi(end,:),[N_x,M]);

subplot(3,2,1)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(2,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
scatter(Xi_est(1,:),Xi_est(2,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$M = 64$');hold on
xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(a) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);

subplot(3,2,2)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(3,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
scatter(Xi_est(1,:),Xi_est(3,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$M = 64$');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(b) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);

load('FPPS_J256.mat');
Xi_est = reshape(Xi(end,:),[N_x,M]);

subplot(3,2,3)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(2,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
scatter(Xi_est(1,:),Xi_est(2,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$M = 128$');hold on
xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(c) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);

subplot(3,2,4)
scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
scatter(Xi_sample(1,:),Xi_sample(3,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
scatter(Xi_est(1,:),Xi_est(3,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$M = 128$');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(d) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);


load('FPPS_J512.mat');
Xi_est = reshape(Xi(end,:),[N_x,M]);

subplot(3,2,5)
s1 = scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(2,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
s2 = scatter(Xi_sample(1,:),Xi_sample(2,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
s3 = scatter(Xi_est(1,:),Xi_est(2,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$M = 256$');hold on
xlabel('x_1')
ylabel('x_2')
xlim([-1.3,-0.5])
ylim([-0.8,0])
str = sprintf('(e) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);

subplot(3,2,6)
s4 = scatter(Xi_MCMC(1,end-20*9999:20:end),Xi_MCMC(3,end-20*9999:20:end),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0.3],'DisplayName','MCMC');hold on
s5 = scatter(Xi_sample(1,:),Xi_sample(3,:),'MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.3 0 0],'DisplayName','resampled');hold on
s6 = scatter(Xi_est(1,:),Xi_est(3,:),'MarkerFaceColor',[0 .7 0],'MarkerEdgeColor',[0 .3 0],'DisplayName','$M = 256$');hold on
xlabel('x_1')
ylabel('x_3')
xlim([-1.3,-0.5])
ylim([0,0.6])
str = sprintf('(f) particle system');
title(str,'FontSize',20)
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',16);

%%% Potential comparison different M

fig3 = figure(3);
clf(fig3)
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);

load('FPPS_J128.mat');

V_t = zeros(1,length([1:10,10:10:length(t)]));
count = 1;
for tk = [1:10,10:10:length(t)]
    Xihelp = reshape(Xi(tk,:),[N_x,M]);
    V = zeros(M,1); 
    for l = 1:M
        h = mvnpdf(Xihelp',Xihelp(:,l)',B)';
        V(l,:) = log(mean(h))+1/2*norm(Gamma^(1/2)\(G(Xihelp(:,l))-y),2)^2+1/2*norm(Gamma0^(1/2)\(Xihelp(:,l)-m0),2)^2;
    end
    V_t(1,count) = mean(V);
    count = count+1;
end

loglog(t([1:10,10:10:length(t)]),V_t,'-','LineWidth',2,'Color',[0.7 0 0],'DisplayName','$M = 128$');hold on
legend('show','interpreter','latex','FontSize',20);

load('FPPS_J256.mat');

V_t = zeros(1,length([1:10,10:10:length(t)]));
count = 1;
for tk = [1:10,10:10:length(t)]
    Xihelp = reshape(Xi(tk,:),[N_x,M]);
    V = zeros(M,1); 
    for l = 1:M
        h = mvnpdf(Xihelp',Xihelp(:,l)',B)';
        V(l,:) = log(mean(h))+1/2*norm(Gamma^(1/2)\(G(Xihelp(:,l))-y),2)^2+1/2*norm(Gamma0^(1/2)\(Xihelp(:,l)-m0),2)^2;
    end
    V_t(1,count) = mean(V);
    count = count+1;
end

loglog(t([1:10,10:10:length(t)]),V_t,'--','LineWidth',2,'Color',[0 0.7 0],'DisplayName','$M = 256$');hold on
legend('show','interpreter','latex','FontSize',20);

load('FPPS_J512.mat');

V_t = zeros(1,length([1:10,10:10:length(t)]));
count = 1;
for tk = [1:10,10:10:length(t)]
    Xihelp = reshape(Xi(tk,:),[N_x,M]);
    V = zeros(M,1); 
    for l = 1:M
        h = mvnpdf(Xihelp',Xihelp(:,l)',B)';
        V(l,:) = log(mean(h))+1/2*norm(Gamma^(1/2)\(G(Xihelp(:,l))-y),2)^2+1/2*norm(Gamma0^(1/2)\(Xihelp(:,l)-m0),2)^2;
    end
    V_t(1,count) = mean(V);
    count = count+1;
end

loglog(t([1:10,10:10:length(t)]),V_t,':','LineWidth',2,'Color',[0 0 0.7],'DisplayName','$M = 512$');hold on
legend('show','interpreter','latex','FontSize',20);

title('(a) potential','FontSize',20)


%%% Parameter estimation

fig4 = figure(4);    
clf(fig4)
set(fig4, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);
load('FPPS_J512.mat');

p3 = plot(xx_b,[0;Uest_MCMC;0],'--','Color',[0 0 0.7],'LineWidth',3,'DisplayName','MCMC');hold on
plot(xx_b(2:end-1),LB_MCMC,'-','color',[0.3 0.3 0.3],'DisplayName','credible set');hold on
plot(xx_b(2:end-1),UB_MCMC,'-','color',[0.3 0.3 0.3],'DisplayName','credible set');hold on

Uest_sample_mean = mean(Uest_sample,2);

SD_sample = sqrt(var(Uest_sample')');
LB_sample = Uest_sample_mean-SD_sample;
UB_sample = Uest_sample_mean+SD_sample;
p2 = plot(xx_b,[0;Uest_sample_mean;0],'-.','Color',[0.7 0 0],'LineWidth',3,'DisplayName','derivative-free');hold on
plot(xx,LB_sample,':','color',[0.3 0.3 0.3],'DisplayName','credible set');hold on
plot(xx,UB_sample,':','color',[0.3 0.3 0.3],'DisplayName','credible set');hold on
fill([xx;xx],[LB_sample';UB_sample'],':','edgecolor',[1 0 0],'LineWidth',0.5);
fill([xx',xx']',[LB_MCMC,UB_MCMC]','-','edgecolor',[0 0 1],'LineWidth',0.5);
p1 = plot(xx_b,[0;utrue;0],'black-','LineWidth',3,'DisplayName','truth');hold on


title('(a) deterministic','FontSize',24)
xlabel('x')
legend('show',[p1,p2,p3],'FontSize',16,'Location','southeast')


