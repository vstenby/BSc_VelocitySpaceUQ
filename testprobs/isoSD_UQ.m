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
%S_noisy = add_noise(S,0.01); 
%Some of the elements in S are negative, and therefore the noise will be
%complex.
S_noisy = S + 0.05*randn(size(S));

S = S_noisy;
A = isoSDA(3,isoSDXinfo,isoSDSinfo);

%% UQ on isoSD case.
alpha = 2.2230e+09; %alpha0 is pretty important and dependent on the noise.

[x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A,S,alpha,100);

%%
showDistribution(std(x_sim,0,2),gridinfo)


