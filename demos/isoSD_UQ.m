%% Bi-Maxwellian UQ Example
% This is a demo of Uncertainty Quantification for the bi-Maxwellian
% fast-ion velocity distribution. 

clear, clc, close all

[A, b, x, L, ginfo] = isoSD();

%Display the true solution
figure
showDistribution(x,ginfo); title('True solution')

%Display nonnegative least-squares solution
figure
showDistribution(TikhNN(A,b,0,[]),ginfo); title('NNLS solution')

%0th order Tikhonov.
disp('Solving 0th order Tikhonov.')
alpha0 = logspace(-10,-5,20);
xalpha0 = TikhNN(A, b, alpha0, []);
[ralpha0, idx0] = relerr(x, xalpha0);
fprintf('Optimal solution: alpha = %.2e, r(alpha) = %.5f\n',alpha0(idx0),ralpha0(idx0))

%1st order Tikhonov.
disp('Solving 1st order Tikhonov.')
alpha1 = logspace(0,10,20); 
xalpha1 = TikhNN(A, b, alpha1, L);
[ralpha1, idx1] = relerr(x, xalpha1);
fprintf('Optimal solution: alpha = %.2e, r(alpha) = %.5f\n',alpha1(idx1),ralpha1(idx1))

nsim = 1000; nburnin = 100;

disp('Running Gibbs Sampler with 0th order prior - this might take a while.')
[xsim0, alphasim0, deltasim0, lambdasim0, info0] = NNHGS(A,b,[],nsim);

disp('Running Gibbs Sampler with 1st order prior - this might take a while.')
[xsim1, alphasim1, deltasim1, lambdasim1, info1] = NNHGS(A,b,L,nsim);

%Convergence plots for 0th order
disp('0th order:')
chain_analysis(deltasim0(nburnin:end),lambdasim0(nburnin:end))

%Convergence plots for 1st order
disp('1st order:')
chain_analysis(deltasim1(nburnin:end),lambdasim1(nburnin:end))

%Plot the regularisation parameters with sample quantiles
%0th order
figure
semilogx(alpha0, ralpha0)
hold on
plot(alpha0(idx0),ralpha0(idx0),'k.', 'MarkerSize', 15) %Optimum
alpha0_mean = mean(alphasim0(nburnin:end)); r0mean = relerr(x,TikhNN(A,b,alpha0_mean));
plot(alpha0_mean, r0mean, 'ro', 'MarkerSize', 10)
qalpha0 = quantile(alphasim0(nburnin:end),[0.025, 0.975]);
errorbar(alpha0_mean, r0mean, ...
         qalpha0(1) - alpha0_mean, ...
         qalpha0(2) - alpha0_mean, ...
         'horizontal','LineWidth',1,'color','r', 'HandleVisibility', 'off')
title('0th order regularisation parameter results')

figure
semilogx(alpha1, ralpha1)
hold on
plot(alpha1(idx1),ralpha1(idx1),'k.', 'MarkerSize', 15) %Optimum
alpha1_mean = mean(alphasim1(nburnin:end)); r1mean = relerr(x,TikhNN(A,b,alpha1_mean,L));
plot(alpha1_mean, r1mean, 'ro', 'MarkerSize', 10)
qalpha1 = quantile(alphasim1(nburnin:end),[0.025, 0.975]);
errorbar(alpha1_mean, r1mean, ...
         qalpha1(1) - alpha1_mean, ...
         qalpha1(2) - alpha1_mean, ...
         'horizontal','LineWidth',1,'color','r', 'HandleVisibility', 'off')
title('1st order regularisation parameter results')

%Plot the credibility bounds
figure
subplot(1,2,1)
showDistribution(mean(xsim0(:,nburnin:end),2),ginfo); title('Sample mean, 0th order prior')
subplot(1,2,2)
showDistribution(cbounds(xsim0(:,nburnin:end)),ginfo); title('95% credibility bounds, 0th order prior')

figure
subplot(1,2,1)
showDistribution(mean(xsim1(:,nburnin:end),2),ginfo); title('Sample mean, 1st order prior')
subplot(1,2,2)
showDistribution(cbounds(xsim1(:,nburnin:end)),ginfo); title('95% credibility bounds, 1st order prior')
