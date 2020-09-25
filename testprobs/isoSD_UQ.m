clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate on this grid
[F, isoSDFinfo] = isoSDF(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, isoSDSinfo] = isoSDS(ustruct,phi);
%S_noisy = add_noise(S,0.01); 
%Some of the elements in S are negative, and therefore the noise will be
%complex.
S_noisy = S + 0.05*randn(size(S));

S = S_noisy;
W = isoSDW(3,isoSDFinfo,isoSDSinfo);

%% UQ on isoSD case.
alpha = 2.2230e+09; %alpha0 is pretty important and dependent on the noise.

[F_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(W,S,alpha,100);

%%
showDistribution(std(F_sim,0,2),gridinfo)


