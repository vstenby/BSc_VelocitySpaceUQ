%biMax Recon. Reconstruction of the drifting bi-Maxwellian distribution.
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid - check construct_vgrid() for default values.
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate the bi-Maxwellian on this grid with default values.
[F, biMaxFinfo] = biMaxF(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, biMaxSinfo] = biMaxS(ustruct,phi);
S_noisy = add_noise(S,0.01);

%Generate A from vpara, vperp, u and phi from biMaxX and biMaxS.
%First argument is ubroadening, which is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
W = biMaxW(3,biMaxFinfo,biMaxSinfo);


%Reconstruct using given alpha-value. 
%This alpha-value is found by trial and error in the biMax_recon.m
%testprob.

%1D:

vpara_linspace = gridinfo.vpara_ax;
vperp_linspace = gridinfo.vperp_ax;
%Another attempt.



%%
%n_x = nvpara;
%n_y = nvperp;

[L1vpara,L1vperp] = gradient_v_space_matrix(vpara_linspace,vperp_linspace,'custom');

x_temp = mosek_TikhNN(W, S, 0);

showDistribution(x_temp,gridinfo)
%L = @(n) spdiags([-ones(n,1) ones(n,1)],[-1 0],n-1,n);
%I = @(n) speye(n,n);
%n = 100;

%n_x = 50;
%n_y = 100;
%L2D = [kron(L(n_x),I(n_y)) ; kron(I(n_x),L(n_y))];

%size(L2D)

%WalphaL=double([transfer_matrix; sqrt(alpha(i))*L1vpara; sqrt(alpha(i))*L1vperp; 2e-6*sqrt(alpha(i))*L0]);
          
%L = [speye(5000) ; L(5000)];
%[~,R] = qr(L);


%L1 = spdiags([-ones(n_x,1) ones(n_x,1)],[-1 0],n_x-1,n_x); 
%I = speye(n_x,n_x); 
%Dx = kron(L1,I);

%D = spdiags([-ones(n_y,1) ones(n_y,1)],[-1 0],n_y+1,n_y); I = speye(n_x,n_x);
%Dy = kron(D,I);

%Lsq = Dx'*Dx + Dy'*Dy;
%L = Lsq' * Lsq;

%alphaN = 100;
%alphavec1 = logspace(-24,-27,alphaN);
%alphavec1 = linspasce(0,1e-27);
%[xalpha, error_dist1] = opt_alpha(alphavec1, W, S_noisy, F(:), Lsq);

%%
semilogx(alphavec1,error_dist1);
xlabel('\alpha'); ylabel('||x_\alpha-x_{true}||^2')

%%
 
%%
%D = spdiags([-ones(nvpara,1) 2*ones(nvpara,1) -ones(nvpara,1) 2*ones(nvpara,1)],[-1:1],nvpara,nvpara);
%L1D_x = spdiags([-ones(nvpara,1) 2*ones(nvpara,1) -ones(nvpara,1) 2*ones(nvpara,1)],[-1:1],nvpara,nvpara);
%L1D_y = spdiags([-ones(nvperp,1) 2*ones(nvperp,1) -ones(nvperp,1) 2*ones(nvperp,1)],[-1:1],nvperp,nvperp);

%L2D_1 = kron(speye(nvpara),L1D_y) + kron(L1D_x,speye(nvperp));
%L2D_2 = kron(L1D_y,speye(nvpara)) + kron(speye(nvperp),L1D_x);

%Lsq       = Ds'*Ds + Dt'*Dt;
%L         = Lsq'*Lsq;        % discrete biharmonic.
    
%%
alphaN = 50;
alphavec1 = logspace(-25,-26,alphaN);
[xalpha, error_dist1] = opt_alpha(alphavec1, W, S_noisy, F(:), Lsq);

semilogx(alphavec1, error_dist1)

%alpha = 2e-26;

%F_recon = mosek_TikhNN(W, S_noisy, alpha, Lsq);

%Display the distribution
%figure(1)
%showDistribution(F, gridinfo)
%title('Bi-Maxwellian distribution')

%figure(2)
%plot(S, 'LineWidth',2)
%hold on
%plot(S_noisy, 'LineWidth',1)
%title('Analytic projection')
%legend('Clean','Noisy')

%figure(3)
%showDistribution(F_recon, gridinfo)
%title('Reconstructed bi-Maxwellian distribution')




