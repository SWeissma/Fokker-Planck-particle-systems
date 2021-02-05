%%% Interacting Langevin particle dynamics without correction

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
End_time = 200; % final time
[Xi_path,time] = Langevin_euler(G,Gamma0,xi0,y,End_time,Gamma);

runningtime = toc;

% computation of spread and residuals of the particle system
spread = zeros(1,length(time));
residual = zeros(1,length(time));
for tk = 1:length(time)
    Xihelp = reshape(Xi_path(tk,:),[I,M2]);
    spread(tk) = 1/M2*norm(Xihelp-mean(Xihelp,2),'fro');
    residual(tk) = 1/M2*norm(Xihelp-Xi_true,'fro');
end

Xi_est = reshape(Xi_path(end,:),[I,M2]);

% collect sample after some iterations
Xi_sample_help = Xi_path(end-10*(2^9/M2)+1:10:end,:);
Xi_sample = reshape(Xi_sample_help(1,:),[I,M2]);
for i = 2:size(Xi_sample_help,1)
    Xi_sample = [Xi_sample, reshape(Xi_sample_help(i,:),[I,M2])];
end

% save the results
save('Langevin_without_correction_M16.mat','runningtime','Xi_sample','spread','residual','time','Xi_est')
% save('Langevin_without_correction_M32.mat','runningtime','Xi_sample','spread','residual','time','Xi_est')
% save('Langevin_without_correction_M64.mat','runningtime','Xi_sample','spread','residual','time','Xi_est')





