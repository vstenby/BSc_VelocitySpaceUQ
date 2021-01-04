%% Bi-Maxwellian Distribution example.
clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[60 80];

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

%%
%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

%Construct the right hand side.
gu = biMaxg(u,phi);
du = u(2)-u(1);
b = du*gu;
[b, e] = add_noise(b);

%Normalize with this noise.
[A, b] = error_normalization(A,b,e);    

%Regularization matrix L (1st order Tikhonov)
L = reguL(vpara, vperp);

%A vector of alphas for finding the optimal regularization parameter.
alpha0 = logspace(-15,-9,50); 
alpha1 = logspace(-4,4,50);
[xalpha0, ~, r0] = TikhNN(A,b,alpha0,[],'return_relerr',true,'x_true',xtrue);
[xalpha1, ~, r1] = TikhNN(A,b,alpha1,L,'return_relerr',true,'x_true',xtrue);

%Find the optimal alpha
[minr0, idx0] = min(r0); alpha0opt = alpha0(idx0);
[minr1, idx1] = min(r1); alpha1opt = alpha1(idx1);

%% Do the Gibbs Sampling
[xNNHGS0, alpha0sim, delta0sim, lambda0sim, info0] = NNHGS(A,b,[],2000,'welfordsalgorithm',true,'nburnin',200);
[xNNHGS1, alpha1sim, delta1sim, lambda1sim, info1] = NNHGS(A,b,L,2000,'welfordsalgorithm',true,'nburnin',200);

%% 0th order bi-Maxwellian
figure
semilogx(alpha0, r0, 'k-')
hold on
%Find the quantiles
CIalpha0 = quantile(alpha0sim,[0.025, 0.975]);
xline(CIalpha0(1),'r-.','LineWidth',3,'alpha',0.2);
xline(mean(alpha0sim),'r-','LineWidth',3,'alpha',0.9);
plot(alpha0opt, minr0, 'b.', 'MarkerSize', 15)
xline(CIalpha0(2),'r-.','LineWidth',3,'alpha',0.2);
xlabel('\alpha','FontSize',25)
ylabel('Relative error','FontSize',15)
legend('r(\alpha)','95% quantiles of \alpha', 'Mean \alpha', 'Optimal \alpha', 'Location','nw','FontSize',15)
title('0th order, bi-Maxwellian')

%% 1st order bi-Maxwellian
figure
semilogx(alpha1, r1, 'k-')
hold on
%Find the quantiles
CIalpha1 = quantile(alpha1sim,[0.025, 0.975]);
xline(CIalpha1(1),'r-.','LineWidth',3,'alpha',0.2);
xline(mean(alpha1sim),'r-','LineWidth',3,'alpha',0.9);
plot(alpha1opt, minr1, 'b.', 'MarkerSize', 15)
xline(CIalpha1(2),'r-.','LineWidth',3,'alpha',0.2);

legend('r(\alpha)','95% quantiles of \alpha', 'Mean \alpha', 'Optimal \alpha', 'Location','nw','FontSize',15)
title('1st order, bi-Maxwellian')
%% Show 0th order
figure
showDistribution(xtrue,ginfo); title('True solution')

figure
showDistribution(xNNHGS0(:,1),ginfo); title('Mean of samples')

figure
showDistribution(xNNHGS0(:,2),ginfo); title('Standard deviation of samples')

figure
showDistribution(TikhNN(A,b,alpha0opt,[]),ginfo); title('Optimal reconstruction')

figure
showDistribution(TikhNN(A,b,mean(alpha0sim)),ginfo); title('Reconstruction using mean alpha')

%% Show 1st order
figure
showDistribution(xtrue,ginfo); title('True solution')

figure
showDistribution(xNNHGS1(:,1),ginfo); title('Mean of samples')

figure
showDistribution(xNNHGS1(:,2),ginfo); title('Standard deviation of samples')

figure
showDistribution(TikhNN(A,b,alpha1opt,L),ginfo); title('Optimal reconstruction')

figure
showDistribution(TikhNN(A,b,mean(alpha1sim),L),ginfo); title('Reconstruction using mean alpha')

