clear, clc, close all

N = 124;
theta = 180-95 : 1.25 : 180+95;         %Angles.
p = 75;                                 %Test setup.
R = 6;
dw = 1.5;
sd = 9;

A = fanlineartomo(N,theta,p,R,dw,sd);

x = dental_phantom(N); 
imagesc(x); axis image

Ax = A*x(:);

SNR = 50;
err_lvl = 100/SNR;

noise = err_lvl/100 * norm(Ax(:)) / sqrt(size(A,1));

b = Ax + noise*randn(size(Ax));

imagesc(reshape(b, p, length(theta)));

%% UQ Using the new tools
L = reguL(124,124);
alphavec = logspace(-1,1,50);
[~,~,r0] = TikhNN(A,b,alphavec,[],'solver','GPCG','return_relerr',true,'x_true',x(:));
[~,~,r1] = TikhNN(A,b,alphavec,L,'solver','GPCG','return_relerr',true,'x_true',x(:));

%% Find the optimal

[~,idx0] = min(r0); optalpha0 = alphavec(idx0);
[~,idx1] = min(r1); optalpha1 = alphavec(idx1);

xopt0 = TikhNN(A,b,optalpha0,[],'solver','GPCG');
xopt1 = TikhNN(A,b,optalpha1,L,'solver','GPCG');

figure
subplot(1,3,1)
imagesc(x); title('True solution'); axis image

subplot(1,3,2)
imagesc(reshape(xopt0,124,124)); title('Optimal 0th order Tikhonov'); axis image;

subplot(1,3,3)
imagesc(reshape(xopt1,124,124)); title('Optimal 1st order Tikhonov'); axis image;

%% Let's go with the UQ

[xNNHGS0, alpha0] = NNHGS(A,b,[],100,'welford',true,'nburnin',10, 'solver', 'GPCG');
[xNNHGS1, alpha1] = NNHGS(A,b,L,100,'welford',true,'nburnin',10, 'solver', 'GPCG');

%%
figure
subplot(1,2,1)
imagesc(reshape(xNNHGS0(:,1),124,124)); axis image;

subplot(1,2,2)
imagesc(reshape(xNNHGS0(:,2),124,124)); axis image; 

figure
subplot(1,2,1)
imagesc(reshape(xNNHGS1(:,1),124,124)); axis image;

subplot(1,2,2)
imagesc(reshape(xNNHGS1(:,2),124,124)); axis image

%%
figure
semilogx(alphavec,r0)
qalpha0 = quantile(alpha0, [0.025, 0.975]);
xline(qalpha0(1),'b')
xline(qalpha0(2),'b')
hold on
semilogx(alphavec,r1)
qalpha1 = quantile(alpha1, [0.025, 0.975]);
xline(qalpha1(1),'b-')
xline(qalpha1(2),'b-')
