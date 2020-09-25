%Example of reconstruction using the noisy analytic signal. 
%In this script, we normalize by the measurement uncertainty and then
%normalize with the 2-norm.

clear, clc, close all

%Add dependencies.
addpath('../functions')
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

%Generate A from vpara, vperp, u and phi from biMaxX and biMaxS.
%First argument is ubroadening, which is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
A = biMaxA(3,biMaxXinfo,biMaxSinfo);

%Normalize by the measurement uncertainty
[Shat, Ahat] = measurement_normalization(S_noisy,A,s);

%Normalize with the 2-norm
[S2hat, A2hat, factor] = numeric_normalization(Shat, Ahat);

%Reconstruct using given alpha-value.
alpha = 3.4e-16;
X_recon = mosek_TikhNN(A2hat, S2hat, alpha);

%Display the distribution
figure(1)
showDistribution(X, gridinfo)
title('Default bi-Maxwellian distribution')

figure(2)
plot(S, 'LineWidth',2)
hold on
plot(S_noisy, 'LineWidth',1)
title('Analytic projection')
legend('Clean','Noisy')

figure(3)
showDistribution(X_recon, gridinfo)
title('Reconstructed bi-Maxwellian distribution')