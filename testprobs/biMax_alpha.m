%Test problem where we find the optimal alpha value by trial and error
%for the bi-Maxwellian case.
clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Default biMax-distribution (x)
[x_true, xinfo] = biMaxx(vpara,vperp);

%Construct the analytic projections (b)
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];
[b, binfo] = biMaxb(ustruct,phi);
[b_noisy, ~, e] = add_noise(b,0.01); 
%Construct the A matrix relating the two.
ubroadening = 3; %default
A = biMaxA(ubroadening,xinfo,binfo);

%Normalization as done in AnalyticTestCase.m
[A_norm, b_norm] = error_normalization(A,b_noisy,e);


%%

%Number of alphas we should be looking at.
alphaN = 20;
alpha  = logspace(-15,-5,alphaN);

%Call mosek_Tikh for all alpha values. This takes a while.
tic
[x, resnorm, residual, exitflag, output] = mosek_TikhNN(A_norm, b_norm, alpha);
toc

[r, min_idx] = relerr(x_true(:),x);

figure
semilogx(alpha,r); title('alpha and relative error, zero-order Tikhonov')
hold on
plot(alpha(min_idx),r(min_idx),'k.','MarkerSize',15)
title(sprintf('alpha=%.2e, relerr = %.4f, zero-order Tikhonov',alpha(min_idx),r(min_idx)))

figure
subplot(1,2,1)
showDistribution(x(:,min_idx),gridinfo); title('xalpha, zero-order Tikhonov')
subplot(1,2,2)
showDistribution(x_true,gridinfo); title('x')


%% Try the different L's. First, we try Per Christians L.
nx = length(gridinfo.vpara_ax);
ny = length(gridinfo.vperp_ax);

%L = reguL(nx,ny);

%ex = ones(nx,1);
%Dxx = spdiags([ex -2*ex ex], [-1 0 1], nx, nx);

%[L1vpara, L1vperp] = gradient_v_space_matrix(gridinfo.vpara_ax, gridinfo.vperp_ax, 'custom');

%Lpch2 = [L1vpara;L1vperp];

%figure
%spy(Lpch)

%figure
%spy(Lpch2)
%LtL = L1vpara'*L1vpara + L1vperp'*L1vperp;
%LtL = sparse(LtL);
%L = chol(LtL);

I = @(n) speye(n);
D = @(n) spdiags([-ones(n,1) ones(n,1)],[-1 0],n,n);
L = [kron(D(ny), speye(nx)) ; kron(speye(ny),D(nx))];
L = chol(L'*L);

%%
%tic
%x0th = mosek_TikhNN(A,b_noisy,1e8);


x1st = mosek_TikhNN(A,b_noisy,1.0e-24,L);
%x0th = mosek_TikhNN(A_norm,b_norm,1e-8);
%x1st = mosek_TikhNN(A_norm,b_norm,1e-8,L);

%figure
%showDistribution(x0th,gridinfo); title('0th order')
figure
showDistribution(x1st,gridinfo); title('1st order')


%%
%BiMax UQ
nsim = 3;
tic
[x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A,b_noisy,1.0e-25,[],nsim,1);
toc

%% Reconstruct using the alpha used.





%analyze_sim(nsim, x_sim, alph_sim, del_sim, lam_sim, gridinfo);

%% Try using this.


%% Alpha values for 1st order Tikhonov.
%Number of alphas we should be looking at.
%alphaN = 5;
%alpha  = logspace(-30,-25,alphaN);

alpha = 1e-45;
L = reguL(gridinfo);

% ~365 seconds when running on my MacBook.
%tic
[x, resnorm, residual, exitflag, output] = mosek_TikhNN(A_norm, b_norm, alpha, L);
%toc

figure
showDistribution(x,gridinfo)

[r, min_idx] = relerr(x_true(:),x);

%% 
[r, min_idx] = relerr(x_true(:),x);

figure
semilogx(alpha,r); title('alpha and relative error, zero-order Tikhonov')
hold on
plot(alpha(min_idx),r(min_idx),'k.','MarkerSize',15)
title(sprintf('alpha=%.2e, relerr = %.4f, 1st order Tikhonov',alpha(min_idx),r(min_idx)))

%Optimal solution
figure
showDistribution(x(:,25),gridinfo)
