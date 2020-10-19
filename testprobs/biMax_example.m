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

example = 2;

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
        
        [vpara, vperp, gridinfo] = construct_vgrid(vparamin,vparamax,vparadim,vperpmin,vperpmax,vperpdim);
        [x_true, xinfo] = biMaxx(vpara,vperp,Tpara,Tperp,vparadrift,options);
        x_true = x_true(:);
        
        [b, binfo] = biMaxb(u,phi,Tpara,Tperp,vparadrift,options);
        
        A = transferMatrix(vpara,vperp,phi,u);
        
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
             
        %Observation angles
        phi=[10 20 40 70 85];
        
        [b, binfo] = biMaxb(ustruct,phi);
        
        u = construct_uvec(ustruct);
        
        A = transferMatrix(vpara,vperp,phi,u);
        
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

%% 2) Reconstruction
x0th = GPCG_TikhNN(A, b_noisy, alphavec);    %0th order Tikhonov
x1st = GPCG_TikhNN(A, b_noisy, alphavec, L); %1st order Tikhonov

%Show the relative error as a function of alpha.
[r0, idx0] = relerr(x_true, x0th);
[r1, idx1] = relerr(x_true, x1st);

semilogx(alphavec,r0, 'r--')
hold on
semilogx(alphavec,r1, 'b--')
hold on

%Plot the mininums
plot(alphavec(idx0),r0(idx0), 'r.', 'MarkerSize',15)
plot(alphavec(idx1),r1(idx1), 'b.', 'MarkerSize',15)
xlabel('\alpha','Fontsize',15)
ylabel('Relative error')
legend('0th order Tikhonov', '1st order Tikhonov')

%Plot the best reconstructions along with the true solution.
figure
showDistribution(x_true,gridinfo); title('True solution');

figure
showDistribution(x0th(:,idx0),gridinfo); title('Optimal 0th order Tikhonov');

figure
showDistribution(x1st(:,idx1),gridinfo); title('Optimal 1st order Tikhonov');

%% 3) Uncertainty Quantification
% 3a) UQ for the 0th order Tikhonov formulation
%[xUQ0th, del0th, lam0th, alph0th] = NNHGS(A, b_noisy, 0, [], nsim);

% 3b) UQ for the 1st order Tikhonov formulation
welford.nburnin = 10;
[xUQ1st, alph1st, del1st, lam1st, info] = NNHGS(A, b_noisy, L, nsim, welford);

%% Analysis
%Analyze simulations
clc
analyze_sim(xUQ1st, alph1st, del1st, lam1st, gridinfo);

%Find the reconstruction with mean alpha.
xalpha = GPCG_TikhNN(A,b_noisy,mean(alph1st),L);


