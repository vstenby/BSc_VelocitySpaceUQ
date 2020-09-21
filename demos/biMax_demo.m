%% Parallel-drifting Bi-Maxwellian example.
clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Set seed for noise.
rng('Default') 

%% Setting up the problem by generating A, X and S. 
% Constructing the (vpara, vperp) grid.
%vparamin, vparamax, vparadim, vperpmin, vperpmax, vperpdim
[vpara, vperp, gridinfo] = construct_vgrid(-1.4e7, 1.4e7, 100, 1e5, 1.4e7, 50);

% Construct our distribution.
[X, biMaxXinfo] = biMaxX(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;

%Number of points per spectrum
ustruct.udim = 200;

%Observation angles
phivec=[10 20 40 70 85];

[S, biMaxSinfo] = biMaxS(ustruct, phivec);

%First input is ubroadening, which is the spectral resolution of the measurements divided by bin width u of the spectra.
A = biMaxA(3, biMaxXinfo, biMaxSinfo);

% Noises is added to the analytic projection.
[S_noisy, s] = add_noise(S,0.01);

%% Find the optimal alpha. 

%This takes a while (~300s on my MacBook)
alphavec1 = logspace(6,9,50);
[~, error_dist1] = opt_alpha(alphavec1, A, S_noisy, X);
%plot(alphavec1, error_dist1)
[~, idx] = min(error_dist1);
alpha_recon1 = alphavec1(idx);

%Optimal alpha I found for a given run with noise_level = 0.01.
%alpha = 2.8118e+08;

X_recon1 = mosek_TikhNN(A, S_noisy, alpha_recon1);



%% Normalizing example.
[Shat, Ahat] = measurement_normalization(S_noisy, A, s);
[S2hat, A2hat, factor] = numeric_normalization(Shat, Ahat);

%This takes a while  (~380s on my Macbook)
alphavec2 = logspace(-17, -14, 50);
[~, error_dist2] = opt_alpha(alphavec2, A2hat, S2hat, X, factor);
%plot(alphavec2, error_dist2)
[~, idx] = min(error_dist2);
alpha_recon2 = alphavec2(idx);

%Optimal alpha I found for a given run.
%alpha = 4.9417e-15;

X_recon2 = mosek_TikhNN(A2hat, S2hat, alpha_recon2).*factor;

%% Show the reconstructions.

figure(1)
showDistribution(X, gridinfo);
title('Original distribution')

figure(2)
showDistribution(X_recon1, gridinfo);
title('Reconstructed distribution without any normalization')

figure(3)
showDistribution(X_recon2, gridinfo);
title('Reconstructed distribution with normalization')

figure(4)
plot(alphavec1, error_dist1)

figure(5)
plot(alphavec2, error_dist2)

%data1 is from non-normalized.
%data2 is from normalized.



