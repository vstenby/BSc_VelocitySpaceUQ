clear, clc, close all

N = 30;
theta = 180-95 : 10 : 180+95;         %Angles.
p = 35;                                 %Test setup.
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

%%
% UQ Using the new tools
L = reguL(N,N);
alphavec = logspace(-2,1,30);
[xrecon,~,r0] = TikhNN(A,b,alphavec,L,'solver','\','return_relerr',true,'x_true',x(:));

semilogx(alphavec,r0)

%%
figure
semilogx(alphavec,r0)

figure
%Optimal recon
[~,idx] = min(r0); imagesc(reshape(xrecon(:,idx),N,N)); axis image;

%% Let's go with the UQ

[xNNHGS0, alpha0] = NNHGS(A,b,[],500,'solver','lsqlin', 'welfordsalgorithm',true,'nburnin',50);

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
