%% Bi-Maxwellian Distribution example.
clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

rhs = 1;

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[10 20 40 70 85];

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;

[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

switch rhs
    case 1
        [b, e] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 
    case 2
        vparadimfine = 100; vperpdimfine = 50;
        [vparafine, vperpfine, ginfofine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
        xfine = biMaxx(vparafine, vperpfine); 
        Afine = transferMatrix(vparafine,vperpfine,phi,u); 
        [b, e] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.
    case 3
        gu = biMaxg(u,phi);
        du = u(2)-u(1);
        b = du*gu;
        [b, e] = add_noise(b);
end

%Normalize with this noise.
[A, b] = error_normalization(A,b,e);    

%Regularization matrix L (1st order Tikhonov)
%L = reguL(vperpdim, vparadim);
L = reguL(vpara,vperp);

%A vector of alphas for finding the optimal regularization parameter.
alphas0 = logspace(-9,-4,50);
alphas1 = logspace(-4,4,50); 

%Find the relative error for 0th order and 1st order.
[~, ~, r0] = TikhNN(A, b, alphas0, [], 'return_relerr', true, 'x_true', xtrue);
[~, ~, r1] = TikhNN(A, b, alphas1,  L, 'return_relerr', true, 'x_true', xtrue);

%Find the alpha with the smallest relative error.
[minr0,idx0] = min(r0); alpha0opt = alphas0(idx0);
[minr1,idx1] = min(r1); alpha1opt = alphas1(idx1);

%Find the corresponding solution.
xopt0 = TikhNN(A,b,alpha0opt);
xopt1 = TikhNN(A,b,alpha1opt,L);

%%
%Do the sampling.
[xNNHGS0, alphasim0, deltasim0, lambdasim0, NNHGS0info] = NNHGS(A,b,[],100,'solver', 'lsqnonneg', 'welfordsalgorithm',true,'nburnin',10);
[xNNHGS1, alphasim1, deltasim1, lambdasim1, NNHGS1info] = NNHGS(A,b,L,100,'solver','lsqnonneg', 'welfordsalgorithm', true, 'nburnin',10);

%Calculate reconstruction with mean(alpha)
xsamplealpha0 = TikhNN(A,b,mean(alphasim0));
xsamplealpha1 = TikhNN(A,b,mean(alphasim1),L);

%Empiric confidence intervals for alpha.
CIalpha0 = quantile(alphasim0,[0.025, 0.975]);
CIalpha1 = quantile(alphasim1,[0.025, 0.975]);

%Calculate the caxis for all solutions and standard deviation of solution.
caxis_mu  = [min([xtrue(:) ; xNNHGS0(:,1) ; xNNHGS1(:,1)]), ...
             max([xtrue(:) ; xNNHGS0(:,1) ; xNNHGS1(:,1)])];
        
caxis_std = [min([xNNHGS0(:,2) ; xNNHGS1(:,2)]), ...
             max([xNNHGS0(:,2) ; xNNHGS1(:,2)])];


figure
subplot(4,3,[1,4,7,10])
semilogx(alphas0,r0, 'r-')
hold on
semilogx(alphas1,r1, 'b-')
ax = gca; ax.YGrid = 'on'; ax.GridLineStyle = '-'; yticks(0:0.025:1);
legend('0th','1st', 'AutoUpdate', 'off');
plot(alpha0opt, minr0,'r.','MarkerSize',15);  xlabel('\alpha','FontSize',15)
plot(alpha1opt, minr1,'b.','MarkerSize',15);

xline(CIalpha0(1),'r-.','LineWidth',3,'alpha',0.2);
xline(CIalpha0(2),'r-.','LineWidth',3,'alpha',0.2);
xline(mean(alphasim0),'r-','LineWidth',3,'alpha',0.9);

xline(CIalpha1(1),'b-.','LineWidth',3,'alpha',0.2);
xline(CIalpha1(2),'b-.','LineWidth',3,'alpha',0.2);
xline(mean(alphasim1),'b-','LineWidth',3,'alpha',0.9);

legend('show')
subplot(4,3,2)
showDistribution(xopt0,ginfo,caxis_mu); title('Optimal alpha (0th)')

subplot(4,3,3)
showDistribution(xopt1,ginfo,caxis_mu); title('Optimal alpha (1st)')

subplot(4,3,5)
showDistribution(xNNHGS0(:,1),ginfo,caxis_mu); title('Posterior mean (0th order)')
subplot(4,3,8)
showDistribution(xsamplealpha0,ginfo,caxis_mu); title('Reconstruction with mean alpha (0th order)')
subplot(4,3,11)
showDistribution(xNNHGS0(:,2),ginfo,caxis_std); title('Standard deviation (0th order)')

subplot(4,3,6)
showDistribution(xNNHGS1(:,1),ginfo,caxis_mu); title('Posterior mean (1st order)')
subplot(4,3,9)
showDistribution(xsamplealpha1,ginfo,caxis_mu); title('Reconstruction with mean alpha (1st order)')
subplot(4,3,12)
showDistribution(xNNHGS1(:,2),ginfo,caxis_std); title('Standard deviation (1st order)')
