%% Drifting bi-Maxwellian UQ
%Commented out so it doesn't reset on HPC.
clear, clc, close all

disp('Starting script')

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate the bi-Maxwellian on this grid with default values.
[X, biMaxXinfo] = biMaxX(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, biMaxSinfo] = biMaxS(ustruct,phi);

[S_noisy, s] = add_noise(S,0.01);

S = S_noisy;

A = biMaxA(3,biMaxXinfo,biMaxSinfo);

%% Normalize the problem.
alpha = 2.8e8;

[x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A,S,alpha,5);
