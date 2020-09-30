%% Try to recreate the "AnalyticTestCase" with my own code.
clear, clc, close all

nsim = 250;

%Here, we use some other values for some reason than the default.
nvpara = 40;
nvperp = 20;
[vpara, vperp, gridinfo] = construct_vgrid(-4e6,4e6,nvpara,1e4,4e6,nvperp);
options.Mi = 2*1.6726e-27; 
[x_true, xinfo] = biMaxF(vpara,vperp,20*1000,20*1000,5e5,options);

%phi = [10 20 40 70 85];
phi = [10 20 40];
u   = [-5:0.1:5]*1e6;

[b, binfo] = biMaxS(u,phi,20*1000,20*1000,5e5,options);    
A = biMaxW(1,xinfo,binfo);

[b_noisy, ~, e] = add_noise(b);

%The normalization used in AnalyticTestCase.m
[A_norm, b_norm] = error_normalization(A, b_noisy, e);

%Per Christians way of constructing L.
L = reguL(nvpara,nvperp);

nalpha = 10;
alpha = logspace(-10,-8,nalpha);

tic
[xalpha, resnorm, residual, exitflag, output] = mosek_TikhNN(A_norm,b_norm,alpha);
toc

% Relative error for each of these alphas.
r = relerr(x_true(:),xalpha);
[ropt, idx] = min(r);
alphaopt = alpha(idx);

semilogx(alpha,r); xlabel('\alpha','FontSize',15); ylabel('Relative error')
hold on
plot(alphaopt,ropt,'k.','MarkerSize',15)

[x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A_norm, b_norm, alphaopt, L, nsim, 1);

%Analyze the simulations.
[x_mu, x_std, p] = analyze_sim(nsim, x_sim, alph_sim, del_sim, lam_sim, gridinfo);


% Relative error
r1 = relerr(x_true(:),x_mu); %Mean of posterior distribution.
r2 = relerr(x_true(:),xalpha(:,idx)); %Optimal xalpha

figure
subplot(1,2,1)
showDistribution(x_mu,gridinfo); title(['Posterior mean, relerr=', num2str(r1)])
subplot(1,2,2)
showDistribution(xalpha(:,idx),gridinfo); title(['optimal xalpha, relerr=', num2str(r2)])
