clear, clc, close all

crime = 0;

%set the grids. dvfine for computation of spectra, dv for inversion
%(vpara,vperp) grid
vparamax=4e6;
vparamin=-4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[10 20 40 70 85];
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the fine forward model.
vparadimfine = 100; vperpdimfine = 50;
[vparafine, vperpfine, ginfofine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
xfine = biMaxx(vparafine, vperpfine); 

%Construct our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

L = reguL(vperpdim, vparadim);

A = transferMatrix(vpara,vperp,phi,u); 
Afine = transferMatrix(vparafine,vperpfine,phi,u); 

if crime == 1
    [b, e] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 
    [A, b] = error_normalization(A,b,e);    %Normalize with this noise.
else
    [b, e] = generate_noisy_b(Afine,xfine); %Basically b = A*x + noise 
    [A, b] = error_normalization(A,b,e);    %Normalize with this noise.
end

alpha_relerr = logspace(-10,1,200);
%Find the relative error for 0th order and 1st order.
[~, ~, r0th] = TikhNN(A, b, alpha_relerr, [], 'relerr', true, 'x_true', xtrue);
[~, ~, r1st] = TikhNN(A, b, alpha_relerr,  L, 'relerr', true, 'x_true', xtrue);

%Calculate the optimal reconstructions for 
[minr0,idx0] = min(r0th); optalpha_0th = alpha(idx0);
[minr1,idx1] = min(r1st); optalpha_1st = alpha(idx1);

[xNNHGS0, alphasim0, deltasim0, lambdasim0, NNHGS0info] = NNHGS(A,b,[],100,'solver', 'lsqnonneg', 'welford',true,'nburnin',250);
[xNNHGS1, alphasim1] = NNHGS(A,b,L,100,'solver','lsqnonneg','welford',true,'nburnin',250);

xopt0 = TikhNN(A,b,optalpha_0th);
xopt1 = TikhNN(A,b,optalpha_1st,L);

xsamplealpha0 = TikhNN(A,b,mean(alphasim0));
xsamplealpha1 = TikhNN(A,b,mean(alphasim1),L);

caxis_mu  = [min([xtrue(:) ; xNNHGS0(:,1) ; xNNHGS1(:,1)]), ...
             max([xtrue(:) ; xNNHGS0(:,1) ; xNNHGS1(:,1)])];
        
caxis_std = [min([xNNHGS0(:,2) ; xNNHGS1(:,2)]), ...
             max([xNNHGS0(:,2) ; xNNHGS1(:,2)])];

figure
subplot(4,3,[1,4,7,10])
semilogx(alpha,r0th, 'r-')
hold on
semilogx(alpha,r1st, 'b-')
legend('0th','1st', 'AutoUpdate', 'off');
plot(optalpha_0th, minr0,'r.','MarkerSize',15);  xlabel('\alpha','FontSize',15)
plot(optalpha_1st, minr1,'b.','MarkerSize',15);
CIalpha0 = quantile(alphasim0,[0.025, 0.975]);

xline(CIalpha0(1),'r-.','LineWidth',3,'alpha',0.2);
xline(CIalpha0(2),'r-.','LineWidth',3,'alpha',0.2);
xline(mean(alphasim0),'r-','LineWidth',3,'alpha',0.9);


CIalpha1 = quantile(alphasim1,[0.025, 0.975]);
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
