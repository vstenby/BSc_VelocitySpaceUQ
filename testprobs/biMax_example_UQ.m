%% Bi-Maxwellian Distribution example.
clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[60 80];

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

gu = biMaxg(u,phi);
du = u(2)-u(1);
b = du*gu;
[b, e] = add_noise(b);

%Normalize with this noise.
[A, b] = error_normalization(A,b,e);    

%Regularization matrix L (1st order Tikhonov)
L = reguL(vperpdim, vparadim);

%A vector of alphas for finding the optimal regularization parameter.
alpha_relerr = logspace(-9,-5,150);

