%biMax Recon. Reconstruction of the drifting bi-Maxwellian distribution.
clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Construct the default grid - check construct_vgrid() for default values.
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
S_noisy = add_noise(S,0.01);

%Generate A from vpara, vperp, u and phi from biMaxX and biMaxS.
%First argument is ubroadening, which is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
W = biMaxW(3,biMaxFinfo,biMaxSinfo);

%Reconstruct using given alpha-value. 
%This alpha-value is found by trial and error in the biMax_recon.m
%testprob.

alpha = 2.8e8;
F_recon = mosek_TikhNN(W, S_noisy, alpha);

%Display the distribution
figure(1)
showDistribution(F, gridinfo)
title('Bi-Maxwellian distribution')

figure(2)
plot(S, 'LineWidth',2)
hold on
plot(S_noisy, 'LineWidth',1)
title('Analytic projection')
legend('Clean','Noisy')

figure(3)
showDistribution(F_recon, gridinfo)
title('Reconstructed bi-Maxwellian distribution')


