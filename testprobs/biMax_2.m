%% BiMax example
%
%  This demo consists of the following three parts:
%
%  1) Setting up the problem.
%
%  2) Reconstruction
%     - with 0th and 1st order Tikhonov
%
%  3) Uncertainty Quantification
%     - with 0th and 1st order Tikhonov

%Clear workspace, console etc.
clear, clc, close all

%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Example 1 is a small scale problem, based on AnalyticTestCase.m
%Example 2 is a bigger scale problem, based on biMaxSlowDown4Viktor.m

example = 1;

%% 1) Setting up the problem
switch example
    case 1
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
        
        %alpha values for 0th order Tikhonov and 1st order Tikhonov.
        alphavec = logspace(5, 15, 50);
        
        %Number of simulations for UQ.
        nsim = 500;
    case 2
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
        
        %alpha values for 0th order Tikhonov and 1st order Tikhonov.
        alphavec = logspace(5, 15, 20);
        
        %Number of simulations for UQ.
        nsim = 100;
    otherwise
        error('Wrong example.')
end

%Add noise to b.
[b_noisy, ~, e] = add_noise(b,0.01);

%This is for the 1st order Tikhonov.
L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

%L has to be square for UQ.
%Tried to outcomment this.
%L = chol(L'*L);


%% NNHGS Test
%welford.nburnin = 300;
%welford.keep = 5;

[x, alpha, delta, lambda, info] = NNHGS(A,b,L,welford);