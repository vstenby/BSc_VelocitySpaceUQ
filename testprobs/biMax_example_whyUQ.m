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

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);
%Regularization matrix L (1st order Tikhonov)
L = reguL(vperpdim, vparadim);

% -- Setup 1 ---
phi=[10 20 30 70 85];
A1 = transferMatrix(vpara,vperp,phi,u);

gu = biMaxg(u,phi);
du = u(2)-u(1);
b = du*gu;
[b, e] = add_noise(b);

%Normalize with this noise.
[A1, b1] = error_normalization(A1,b,e);    

% -- Setup 2 ---
phi=[10 20 40 70 85];
A2 = transferMatrix(vpara,vperp,phi,u);

gu = biMaxg(u,phi);
du = u(2)-u(1);
b = du*gu;
[b, e] = add_noise(b);

%Normalize with this noise.
[A2, b2] = error_normalization(A2,b,e); 

% Sample the two posteriors
xsetup1 = NNHGS(A1,b1,L,350,'welfordsalgorithm',true,'nburnin',35);
xsetup2 = NNHGS(A2,b2,L,350,'welfordsalgorithm',true,'nburnin',35);

%% 

figure
subplot(1,2,1)
showDistribution(xsetup1(:,1),ginfo)
subplot(1,2,2)
showDistribution(xsetup2(:,1),ginfo)
%Basically the same

caxis_std = [min([xsetup1(:,2) ; xsetup2(:,2)]); ...
             max([xsetup1(:,2) ; xsetup2(:,2)])];
figure
subplot(1,2,1)
showDistribution(xsetup1(:,2),ginfo,caxis_std)
subplot(1,2,2)
showDistribution(xsetup2(:,2),ginfo,caxis_std)

