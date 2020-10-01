%% Try to recreate the "AnalyticTestCase" with my own code.
clear, clc, close all

nsim = 1000;

%Here, we use some other values for some reason than the default.
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

%Alphas
nalpha = 10;
alpha = logspace(-10,-5,nalpha);

[vpara, vperp, gridinfo] = construct_vgrid(vparamin,vparamax,vparadim,vperpmin,vperpmax,vperpdim);
[x_true, xinfo] = biMaxx(vpara,vperp,Tpara,Tperp,vparadrift,options);
[b, binfo] = biMaxb(u,phi,Tpara,Tperp,vparadrift,options);    
A = biMaxA(ubroadening,xinfo,binfo);

[b_noisy, ~, e] = add_noise(b); %Add noise to b.

%The normalization used in AnalyticTestCase.m
[A_norm, b_norm] = error_normalization(A, b_noisy, e);

%Per Christians way of constructing L.
L = reguL(vparadim,vperpdim);


%%
[L1vpara, L1vperp] = gradient_v_space_matrix(gridinfo.vpara_ax, gridinfo.vperp_ax, 'custom');
LtL = L1vpara'*L1vpara + L1vperp'*L1vperp;




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
