%Test problem where we find the optimal alpha value by trial and error
%for the bi-Maxwellian case.
clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Run biMax_recon2 example to get the different variables.
run('biMax_recon2.m'); close all
clearvars -except A A2hat S_noisy S2hat factor X gridinfo
fprintf('Variables loaded from reconstruction example 2.\n\n')

alphaN = 50; %Number of alphas to be checked in the linspace.
%50 takes approx. 11 minutes (on my MacBook)

%In both cases, we need X to compare our reconstruction.
%To find the optimal alpha for the non-normalized case, we need:
%A, S_noisy

%To find the optimal alpha for the normalized case, we need:
%A2hat (normalized A), S2hat (normalized S), factor (to get to the absolute
%valued case)

% Find the optimal alpha for the normal case.
%This takes a while (~300s on my MacBook)
fprintf('Finding optimal alpha for the normal case.\n')
alphavec1 = logspace(6,9,alphaN);
[~, error_dist1] = opt_alpha(alphavec1, A, S_noisy, X);
%plot(alphavec1, error_dist1)
[~, idx] = min(error_dist1);
alpha_recon1 = alphavec1(idx);
fprintf('Optimal alpha found.\n\n')

%Optimal alpha I found for a given run with noise_level = 0.01.
%alpha = 2.8118e+08;

X_recon1 = mosek_TikhNN(A, S_noisy, alpha_recon1);

% Normalized example.
%This takes a while  (~380s on my Macbook)
fprintf('Finding optimal alpha for the normalized case.\n')
alphavec2 = logspace(-17, -14, alphaN);
[~, error_dist2] = opt_alpha(alphavec2, A2hat, S2hat, X, factor);
%plot(alphavec2, error_dist2)
[~, idx] = min(error_dist2);
alpha_recon2 = alphavec2(idx);
fprintf('Optimal alpha found.\n\n')

%Optimal alpha I found for a given run.
%alpha = 4.9417e-15;

X_recon2 = mosek_TikhNN(A2hat, S2hat, alpha_recon2).*factor;

figure(1)
sgtitle('Finding $\alpha$ for case with and without normalization','interpreter','latex')
subplot(1,2,1)
semilogx(alphavec1, error_dist1)
title('Without normalization'); xlabel('\alpha');
subplot(1,2,2)
semilogx(alphavec2, error_dist2); xlabel('\alpha')
title('With normalization');

