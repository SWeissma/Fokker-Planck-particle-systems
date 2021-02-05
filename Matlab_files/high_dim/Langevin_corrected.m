%%% Interacting Langevin particle dynamics

clear;
load('reference_data.mat')

% path to functions
path('./Functions',path);

% Parameters
M2 = 2^4;   % ensemble size   
% M2 = 2^5;
% M2 = 2^6;

% prior sample
xi0 = mvnrnd(m0,Gamma0,M2)';  % initial ensemble

tic;
End_time = 200;  % number of iterations
[Xi_path2,time2] = Langevin_correction_euler(G,Gamma0,xi0,y,End_time,Gamma);

runningtime2 = toc;

% computation of spread and residuals of the particle system
spread2 = zeros(1,length(time2));
residual2 = zeros(1,length(time2));
for tk = 1:length(time2)
    Xihelp = reshape(Xi_path2(tk,:),[I,M2]);
    spread2(tk) = 1/M2*norm(Xihelp-mean(Xihelp,2),'fro')^2;
    residual2(tk) = 1/M2*norm(Xihelp-Xi_true,'fro')^2;
end

Xi_est2 = reshape(Xi_path2(end,:),[I,M2]);

% collect sample after some iterations
Xi_sample_help = Xi_path2(end-10*(2^9/M2)+1:10:end,:);
Xi_sample2 = Xi_est2;
for i = 2:size(Xi_sample_help,1)
    Xi_sample2 = [Xi_sample2, reshape(Xi_sample_help(i,:),[I,M2])];
end

% save the results
save('Langevin_M16.mat','runningtime2','Xi_sample2','spread2','residual2','time2','Xi_est2')
% save('Langevin_M32.mat','runningtime2','Xi_sample2','spread2','residual2','time2','Xi_est2')
% save('Langevin_M64.mat','runningtime2','Xi_sample2','spread2','residual2','time2','Xi_est2')

