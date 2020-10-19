%% biMax Nsim Simulation 14-10-2020
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

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

%L regularization matrix
L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

%Observation angles
phi = [10 20 40 70 85];

for i=4:5
    nsim = 1e4-1000;
    phitemp = phi(randperm(i));
    disp(phitemp)
    tic
    [b, binfo] = biMaxb(ustruct,phi);
    [b_noisy, ~, e] = add_noise(b,0.01);
    A = biMaxA(ubroadening, xinfo, binfo);
    welford.keep = 1;
    welford.nburnin = 1000;
    [x, alpha, delta, lambda, info] = NNHGS(A,b_noisy,L,nsim, welford);
    save(strcat('./nsim/nangle',num2str(i),'.mat'),'x', 'alpha','delta','lambda','info','phitemp');
    toc
end
