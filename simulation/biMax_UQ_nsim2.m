%% Simulation script for testing the importance of nsim.
clear, clc, close all

%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%% Setting up the problem.
%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();
vparadim = length(gridinfo.vpara_ax);
vperpdim = length(gridinfo.vperp_ax);

%Evaluate the bi-Maxwellian on this grid with default values.
[x_true, xinfo] = biMaxx(vpara, vperp); x_true = x_true(:);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;

%ubroadening is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
ubroadening = 3; 

%Observation angles
phi=[10 20 40 70 85];

[b, binfo] = biMaxb(ustruct,phi);

A = biMaxA(ubroadening, xinfo, binfo);

%Add noise to b.
[b_noisy, ~, e] = add_noise(b,0.01);

%This is for the 1st order Tikhonov.
L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

options.welford = 1;
options.waitbar = 1;
    
%% Setting nsim and folder

folderName = 'nsimUQ2';
mkdir(folderName)

expo = 2:6;
nsim = logspace(expo(1),expo(end),5);

for i=1:length(nsim)
    n = nsim(i);
    simname = strcat(folderName,'/n_1e',num2str(expo(i)),'.mat');
    %Set the same seed
    rng('default')
    [x, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b_noisy, 0, L, n, options);
    save(simname,'x','del_sim','lam_sim','alph_sim', 'n', 'gridinfo')
end

