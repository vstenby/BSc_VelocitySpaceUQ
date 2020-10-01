%% biMax_recon : Experimental Setup 1
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid - check construct_vgrid() for default values.
[vpara, vperp, gridinfo] = construct_vgrid();
vparadim = length(gridinfo.vpara_ax);
vperpdim = length(gridinfo.vperp_ax);

%Evaluate the bi-Maxwellian on this grid with default values.
[x_true, xinfo] = biMaxx(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[b, binfo] = biMaxb(ustruct,phi);
[b_noisy, ~, e] = add_noise(b,0.01);

%Generate A from vpara, vperp, u and phi from biMaxx and biMaxb.
%First argument is ubroadening, which is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
A = biMaxA(3,xinfo,binfo);

[A_norm, b_norm] = error_normalization(A,b_noisy,e);

L = reguL(vparadim,vperpdim);

%Regularization L has to be square for UQ.
L = chol(L'*L);

%% Find xalpha for different alpha values.

alpha_0th = logspace(-10, -8, 20);
alpha_1st = logspace(-10, -8, 20);

%Found by testing
disp('Finding solutions for 0th order Tikhonov.')
x0th = mosek_TikhNN(A_norm,b_norm,alpha_0th);

disp('Finding solutions for 1st order Tikhonov.')
x1st = mosek_TikhNN(A_norm,b_norm,alpha_1st,L);

[r_0th, idx_0th] = relerr(x_true(:),x0th);
[r_1st, idx_1st] = relerr(x_true(:),x1st);

figure
semilogx(alpha_0th, r_0th);
hold on
semilogx(alpha_1st, r_1st);
hold on
semilogx(alpha_0th(idx_0th), r_0th(idx_0th), 'k.','MarkerSize',15)
hold on
semilogx(alpha_1st(idx_1st), r_1st(idx_1st), 'kx','MarkerSize',15)
legend('0th order', '1st order', 'Minimum alpha for 0th', 'Minimum alpha for 1st')

figure
showDistribution(x0th(:,idx_0th),gridinfo); title('0th order Tikhonov solution');

figure
showDistribution(x1st(:,idx_1st),gridinfo); title('1st order Tikhonov solution');

figure
showDistribution(x_true,gridinfo); title('True solution')

%% The best alpha values
alphaopt_0th = alpha_0th(idx_0th); 
alphaopt_1st = alpha_1st(idx_1st); 
%Try and solve using GPCG and the same approach. 
LtL = L'*L;

%% Try and solve the same problem 2 ways.
[~,N] = size(A);
alphas = logspace(5,15,20);
x0th = GPCG_TikhNN(A,b_noisy,alphas);
[r0th, idx0] = relerr(x_true(:),x0th);
r0opt = r0th(idx0);

x1st = GPCG_TikhNN(A,b_noisy,alphas,L);
[r1st, idx1] = relerr(x_true(:),x1st);
r1opt = r1st(idx1);

%% Plot the optimal
figure
showDistribution(x0th(:,idx0),gridinfo); title('0th order Tikhonov')

figure
showDistribution(x1st(:,idx1),gridinfo); title('1st order Tikhonov')

%% UQ with 0th order Tikhonov assumptions.
n_sim = 100;
[x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b_noisy, 0, [], n_sim, 1);

[~, ~, ~, q_alpha] = analyze_sim(100, x_sim, alph_sim, del_sim, lam_sim, gridinfo);

semilogx(alphas,r0th)
hold on
xline(q_alpha(1), 'r--')
xline(q_alpha(2), 'r')
xline(q_alpha(3), 'k')
xline(q_alpha(4), 'b')
xline(q_alpha(5), 'b--')

%% UQ with 1st order Tikhonov assumptions.
n_sim = 100;

[x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b_noisy, 0, L, n_sim, 1);

[~,~,~,q_alpha] = analyze_sim(100, x_sim, alph_sim, del_sim, lam_sim, gridinfo);

%Plot the points
semilogx(alphas,r1st)
hold on
xline(q_alpha(1), 'r--')
xline(q_alpha(2), 'r')
xline(q_alpha(3), 'k')
xline(q_alpha(4), 'b')
xline(q_alpha(5), 'b--')