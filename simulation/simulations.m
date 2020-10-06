%% In this MATLAB script, the different simulation folders can be looked at.
clear, clc, %close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

folder = '10  20  70';
% Analyze the simulation
load(strcat('biMaxUQ_angles/', folder, '/setup.mat'))
load(strcat('biMaxUQ_angles/', folder, '/0thTikhUQ.mat'))
load(strcat('biMaxUQ_angles/', folder, '/1stTikhUQ.mat'))

input('Analyze 0th order Tikhonov ')
disp(' ')
analyze_sim(nsim, xUQ0th, alph0th, del0th, lam0th, gridinfo);
disp(' ')
input('Analyse 1st order Tikhonov ')
analyze_sim(nsim, xUQ1st, alph1st, del1st, lam1st, gridinfo);

