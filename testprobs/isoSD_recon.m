clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate on this grid
[X, isoSDXinfo] = isoSDX(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, isoSDSinfo] = isoSDS(ustruct,phi);

%Because some of the elements in S are negative, this will add complex
%noise. I will have to ask Mirko about this.
%S_noisy = add_noise(S,0.01);
S_noisy = S + 0.05*randn(size(S));

A = isoSDA(3,isoSDXinfo,isoSDSinfo);

%%
%Reconstruct using given alpha-value.
alpha = 2.2230e+09;
X_recon = mosek_TikhNN(A, S_noisy, alpha);

figure
subplot(1,2,1)
showDistribution(X, gridinfo)

subplot(1,2,2)
showDistribution(X_recon,gridinfo)
