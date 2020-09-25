%Test problem where we find the optimal alpha value by trial and error
%for the bi-Maxwellian case.
clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Run isoSD_recon to get the variables.
run('isoSD_recon.m'); close all
clearvars -except A X S_noisy gridinfo
fprintf('Variables loaded from reconstruction example 1.\n\n')
tic
alphaN = 50; %Number of alphas to be checked in the linspace.

% Find the optimal alpha for the normal case.
%This takes a while (~300s on my MacBook)
fprintf('Finding optimal alpha for the normal case.\n')
alphavec1 = logspace(6,10,alphaN);
[~, error_dist1] = opt_alpha(alphavec1, A, S_noisy, X);
plot(alphavec1, error_dist1)
[~, idx] = min(error_dist1);
alpha_recon1 = alphavec1(idx);
fprintf('Optimal alpha found.\n\n')
toc

%%
%Optimal alpha I found for a given run with noise_level = 0.01.
alpha = 2.6827e+09;

X_recon1 = mosek_TikhNN(A, S_noisy, alpha_recon1);



