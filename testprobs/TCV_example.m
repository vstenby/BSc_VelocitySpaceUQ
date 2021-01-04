%% Initialisation
clear, clc, close all

% -- Load the data -- 
load('../data/TCV/ts4.mat');
A = transf_matrix; b = double(S_TR_invcrime_n); xtrue = double(f_TR_coarse(:));

%clearvars -except A b xtrue

%% 
clc

L = reguL(20,20);

alpha = logspace(-20,-10,50);
[x, ~, r] = TikhNN(A,b,alpha,L,'return_relerr', true, 'x_true', xtrue);

%%
figure
subplot(1,2,1)
semilogx(alpha, r)
subplot(1,2,2)
[~,idx] = min(r);
showDistribution(x(:,idx),[20,20])

%%
clear, clc, close all

% -- Load the data -- 
load('../data/TCV/ts4.mat');
A = transf_matrix_broad; b = double(S_data); err = err_data; 
xtrue = double(f_TR_coarse(:));

[A, b] = error_normalization(A,b,err);

clearvars -except A b xtrue 

L = reguL(20,20);

alpha = logspace(-5,5,50);
%[xalpha, ~, r] = TikhNN(A,b,alpha,L,'scaling','norm', 'return_relerr', true, 'x_true', xtrue);

nburnin = 100;
nsim = 1000;
[xsim, alphasim] = NNHGS(A,b,L,nsim,'scaling','norm','welfordsalgorithm',true,'nburnin',nburnin);

%%
figure
semilogx(alpha,r)
qalpha = quantile(alphasim,[0.025, 0.975]);
xline(qalpha(1),'k--')
xline(mean(alphasim),'k-')
xline(qalpha(2),'k--')

figure
subplot(1,2,1)
showDistribution(xtrue,[20,20]);
subplot(1,2,2)
showDistribution(xsim(:,1),[20,20]);



