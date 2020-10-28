%% Third time's the charm...
clear, clc, close all

%set the grids. dvfine for computation of spectra, dv for inversion
%(vpara,vperp) grid
vparamax=4e6;
vparamin=-4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[10 20 40 70 85];

%Add the dependencies
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

[vparafine, vperpfine, ginfofine] = construct_vgrid(100,50,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
[vpara, vperp, ginfo] = construct_vgrid(40,20,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);

xtrue = biMaxx(vpara,vperp);
xfine = biMaxx(vparafine, vperpfine); %Basically the same now.
showDistribution(xfine,ginfofine)

u = construct_uvec('umin', [], 'umax', 5*1e6, 'du', 1e5); %Basically the same.

Afine = transferMatrix(vparafine,vperpfine,phi,u); %1/du * Afine and you get the transfermatrixCTSfine.
A = transferMatrix(vpara,vperp,phi,u); %1/du * A and you get tranfermatrixCTS.

[b, e] = generate_noisy_b(Afine,xfine); 
[A, b] = error_normalization(A,b,e);

alpha = logspace(-10,1,200);
L = reguL(20,40);

[~, ~, r0th] = TikhNN(A, b, alpha, [], 'relerr', true, 'x_true', xtrue);
[~, ~, r1st] = TikhNN(A, b, alpha,  L, 'relerr', true, 'x_true', xtrue);

[xNNHGS0, alphasim0] = NNHGS(A,b,[],100,'solver','GPCG','welford',true,'nburnin',250);
[xNNHGS1, alphasim1] = NNHGS(A,b,L,100,'solver','GPCG','welford',true,'nburnin',250);

figure
subplot(4,3,[1,4,7,10])
semilogx(alpha,r0th, 'r-')
hold on
semilogx(alpha,r1st, 'b-')
legend('0th','1st', 'AutoUpdate', 'off');


[minr0,idx0] = min(r0th); optalpha_0th = alpha(idx0);
[minr1,idx1] = min(r1st); optalpha_1st = alpha(idx1);

hold on
plot(optalpha_0th, minr0,'r.','MarkerSize',15);  xlabel('\alpha','FontSize',15)
plot(optalpha_1st, minr1,'b.','MarkerSize',15);
qalpha0 = quantile(alphasim0,[0 0.25 0.5 0.75 1]);
xline(qalpha0(1),'r-.','LineWidth',3,'alpha',0.2);
xline(qalpha0(2),'r-','LineWidth',3,'alpha',0.2);
xline(qalpha0(3),'r-','LineWidth',3,'alpha',0.9);
xline(qalpha0(4),'r-','LineWidth',3,'alpha',0.2);
xline(qalpha0(5),'r-.','LineWidth',3,'alpha',0.2);

qalpha1 = quantile(alphasim1,[0 0.25 0.5 0.75 1]);
xline(qalpha1(1),'b-.','LineWidth',3,'alpha',0.2);
xline(qalpha1(2),'b-','LineWidth',3,'alpha',0.2);
xline(qalpha1(3),'b-','LineWidth',3,'alpha',0.9);
xline(qalpha1(4),'b-','LineWidth',3,'alpha',0.2);
xline(qalpha1(5),'b-.','LineWidth',3,'alpha',0.2);


legend('show')
subplot(4,3,[2 3])
showDistribution(xtrue,ginfo); title('True solution')

subplot(4,3,5)
showDistribution(xNNHGS0(:,1),ginfo); title('Posterior mean (0th order)')
subplot(4,3,8)
xsamplealpha0 = TikhNN(A,b,optalpha_0th);
showDistribution(xsamplealpha0,ginfo); title('Reconstruction with mean alpha (0th order)')
subplot(4,3,11)
showDistribution(xNNHGS0(:,2),ginfo); title('Standard deviation (0th order)')

subplot(4,3,6)
showDistribution(xNNHGS1(:,1),ginfo); title('Posterior mean (1st order)')
subplot(4,3,9)
showDistribution(xNNHGS1(:,1),ginfo); title('Reconstruction with mean alpha (1st order)')
subplot(4,3,12)
showDistribution(xNNHGS1(:,2),ginfo); title('Standard deviation (1st order)')
