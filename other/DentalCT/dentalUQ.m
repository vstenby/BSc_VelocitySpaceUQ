clear, clc, close all

N = 30;
theta = 180-95 : 10 : 180+95;         %Angles.
p = 35;                               %Test setup.
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
I = @(d) speye(d);
D = @(d) spdiags([zeros(d-1,1) -ones(d-1,1) ones(d-1,1)],-1:1,d-1,d);

L = [kron(D(30), I(30)) ; kron(I(30),D(30))];

alphavec = logspace(-2,1,30);
[x0,~,r0] = TikhNN(A,b,alphavec,[],'return_relerr',true,'x_true',x(:));
[x1,~,r1] = TikhNN(A,b,alphavec,L,'return_relerr',true,'x_true',x(:));

%% Tikhonov

%Relative error
figure
semilogx(alphavec,r0);
[r0min, idx0] = min(r0); hold on; plot(alphavec(idx0),r0min,'k.','MarkerSize',15)
hold on
semilogx(alphavec,r1);
[r1min, idx1] = min(r1); hold on; plot(alphavec(idx1),r1min,'k.','MarkerSize',15)

%Optimal reconstructions
figure
subplot(1,2,1)
imagesc(reshape(x0(:,idx0),N,N)); axis image; title(sprintf('Optimal reconstruction, 0th order. r = %.2f',r0min))

subplot(1,2,2)
imagesc(reshape(x1(:,idx1),N,N)); axis image; title(sprintf('Optimal reconstruction, 1st order. r = %.2f',r1min))

%% Uncertainty Quantification

[xNNHGS0, alpha0] = NNHGS(A,b,[],500,'welfordsalgorithm',true,'nburnin',50);
[xNNHGS1, alpha1] = NNHGS(A,b,L,500,'welfordsalgorithm',true,'nburnin',50);

%%
figure
imagesc(reshape(xNNHGS0(:,1),30,30)); axis image; colorbar()

figure
imagesc(reshape(xNNHGS0(:,2),30,30)); axis image; colorbar()

figure
imagesc(reshape(xNNHGS1(:,1),30,30)); axis image; colorbar()

figure
imagesc(reshape(xNNHGS1(:,2),30,30)); axis image; colorbar()

%%
figure
xalpha_mean = TikhNN(A,b,mean(alpha0));
semilogx(alphavec,r0, 'r-')
hold on
qalpha0 = quantile(alpha0, [0.025, 0.975]);
plot(mean(alpha0), relerr(x, xalpha_mean), 'ro', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alpha0),relerr(x, xalpha_mean), ...
         qalpha0(1)-mean(alpha0), qalpha0(2)-mean(alpha0), ...
         'horizontal','LineWidth',1,'color','r', 'HandleVisibility', 'off')
xlabel('\alpha','FontSize',20); ylabel('Relative error','FontSize',20)
%Plot the found minimum
plot(alphavec(idx0),r0min,'k.','MarkerSize',15, 'HandleVisibility', 'off')

semilogx(alphavec,r1, 'b-')
hold on
qalpha1 = quantile(alpha1, [0.025, 0.975]);
plot(mean(alpha1), relerr(x, xalpha_mean1), 'bo', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alpha1),relerr(x, xalpha_mean1), ...
         qalpha1(1)-mean(alpha1), qalpha1(2)-mean(alpha1), ...
         'horizontal','LineWidth',1,'color','b', 'HandleVisibility', 'off')
xlabel('\alpha','FontSize',20); ylabel('Relative error','FontSize',20)
plot(alphavec(idx1),r1min,'k.','MarkerSize',15)
xlim([1e-2, 1e0])
legend('0th order prior', '1st order prior','Location','nw','FontSize',20)