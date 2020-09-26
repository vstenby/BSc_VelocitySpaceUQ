%% Drifting bi-Maxwellian UQ
%Commented out so it doesn't reset on HPC.
clear, clc, close all

%Set the simulation setup.
nsim = 20; 
disp_waitbar = 1;

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate the bi-Maxwellian on this grid with default values.
[F, biMaxFinfo] = biMaxF(vpara, vperp);

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

W = biMaxW(3,biMaxFinfo,biMaxSinfo);

alpha = 2.8e8; %alpha0

[F_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(W,S,alpha,nsim,disp_waitbar);

analyze_sim(nsim, F_sim, alph_sim, del_sim, lam_sim, gridinfo);

