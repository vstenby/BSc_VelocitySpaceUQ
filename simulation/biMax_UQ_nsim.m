%% Simulation script for testing the importance of nsim.
clear, clc, close all

%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%% Setting up the problem.
vparamin = -4e6;
vparamax = 4e6;
vparadim = 40;

vperpmin = 1e4;
vperpmax = 4e6;
vperpdim = 20;

Tpara = 2e4;
Tperp = 2e4;
vparadrift = 5e5;

options.Mi = 2*1.6726e-27;

phi = [10 20 40 70 85];
u   = [-5e6:1e5:5e6];
ubroadening = 1;

[vpara, vperp, gridinfo] = construct_vgrid(vparamin,vparamax,vparadim,vperpmin,vperpmax,vperpdim);
[x_true, xinfo] = biMaxx(vpara,vperp,Tpara,Tperp,vparadrift,options);
x_true = x_true(:);

[b, binfo] = biMaxb(u,phi,Tpara,Tperp,vparadrift,options);

A = biMaxA(ubroadening, xinfo, binfo);

%Add noise to b.
[b_noisy, ~, e] = add_noise(b,0.01);

%This is for the 1st order Tikhonov.
L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

options.welford = 1;
options.waitbar = 1;
    
%% Setting nsim and folder

folderName = 'nsimUQ';
mkdir(folderName)

expo = 2:7;
nsim = logspace(expo(1),expo(end),6);

for i=1:length(nsim)
    n = nsim(i);
    simname = strcat(folderName,'/n_1e',num2str(expo(i)),'.mat');
    [x, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b_noisy, 0, L, n, options);
    save(simname,'x','del_sim','lam_sim','alph_sim', 'n', 'gridinfo')
end



