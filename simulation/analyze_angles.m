%% Script made for analyzing the angles.
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Load one experimental setup (10,20,40)
load('biMaxUQ_angles/10  20  40/setup.mat')
load('biMaxUQ_angles/10  20  40/0thTikhUQ.mat')
load('biMaxUQ_angles/10  20  40/1stTikhUQ.mat')

alph0th_setup1 = alph0th;
alph1st_setup1 = alph1st;
del0th_setup1 = del0th;
del1st_setup1 = del1st;
lam0th_setup1 = lam0th;
lam1st_setup1 = lam1st;
xUQ0th_setup1 = xUQ0th;
xUQ1st_setup1 = xUQ1st;

%Load the other experimental setup (10,20,70)
load('biMaxUQ_angles/10  20  70/setup.mat')
load('biMaxUQ_angles/10  20  70/0thTikhUQ.mat')
load('biMaxUQ_angles/10  20  70/1stTikhUQ.mat')

alph0th_setup2 = alph0th;
alph1st_setup2 = alph1st;
del0th_setup2 = del0th;
del1st_setup2 = del1st;
lam0th_setup2 = lam0th;
lam1st_setup2 = lam1st;
xUQ0th_setup2 = xUQ0th;
xUQ1st_setup2 = xUQ1st;

clearvars -except alph0th_setup1 alph0th_setup2 ...
                  alph1st_setup1 alph1st_setup2 ...
                  del0th_setup1 del0th_setup2 ...
                  del1st_setup1 del1st_setup2 ...
                  lam0th_setup1 lam0th_setup2 ...
                  lam1st_setup1 lam1st_setup2 ...
                  xUQ0th_setup1 xUQ0th_setup2 ...
                  xUQ1st_setup1 xUQ1st_setup2 ...
                  xtrue gridinfo nsim x_true
              
[xmu_0th_setup1, xstd_0th_setup1] = analyze_sim(nsim, xUQ0th_setup1, alph0th_setup1, del0th_setup1, lam0th_setup1, gridinfo);
[xmu_1st_setup1, xstd_1st_setup1] = analyze_sim(nsim, xUQ1st_setup1, alph1st_setup1, del1st_setup1, lam1st_setup1, gridinfo);                  
[xmu_0th_setup2, xstd_0th_setup2] = analyze_sim(nsim, xUQ0th_setup2, alph0th_setup2, del0th_setup2, lam0th_setup2, gridinfo);
[xmu_1st_setup2, xstd_1st_setup2] = analyze_sim(nsim, xUQ1st_setup2, alph1st_setup2, del1st_setup2, lam1st_setup2, gridinfo);

%% Show the four means with the same colorscale
close all
mucaxis(1) = min([xmu_0th_setup1 ; xmu_0th_setup2 ; xmu_1st_setup1 ; xmu_1st_setup2]);
mucaxis(2) = max([xmu_0th_setup1 ; xmu_0th_setup2 ; xmu_1st_setup1 ; xmu_1st_setup2]);

figure
sgtitle('Mean of 18000 samples from posterior')
subplot(2,2,1)
showDistribution(xmu_0th_setup1,gridinfo); caxis(mucaxis); title('phi = [10 20 40], 0th order nn. Tikh')
subplot(2,2,2)
showDistribution(xmu_0th_setup2,gridinfo); caxis(mucaxis); title('phi = [10 20 70], 0th order nn. Tihk')
subplot(2,2,3)
showDistribution(xmu_1st_setup1,gridinfo); caxis(mucaxis); title('phi = [10 20 40], 1st order nn. Tikh')
subplot(2,2,4)
showDistribution(xmu_1st_setup2,gridinfo); caxis(mucaxis); title('phi = [10 20 70], 1st order nn. Tikh')
       
[xmu_0th_setup1, xstd_0th_setup1] = analyze_sim(nsim, xUQ0th_setup1, alph0th_setup1, del0th_setup1, lam0th_setup1, gridinfo);
[xmu_1st_setup1, xstd_1st_setup1] = analyze_sim(nsim, xUQ1st_setup1, alph1st_setup1, del1st_setup1, lam1st_setup1, gridinfo);                  
[xmu_0th_setup2, xstd_0th_setup2] = analyze_sim(nsim, xUQ0th_setup2, alph0th_setup2, del0th_setup2, lam0th_setup2, gridinfo);
[xmu_1st_setup2, xstd_1st_setup2] = analyze_sim(nsim, xUQ1st_setup2, alph1st_setup2, del1st_setup2, lam1st_setup2, gridinfo);

%% Show the four means with the same colorscale
close all
mucaxis(1) = min([xmu_0th_setup1 ; xmu_0th_setup2 ; xmu_1st_setup1 ; xmu_1st_setup2]);
mucaxis(2) = max([xmu_0th_setup1 ; xmu_0th_setup2 ; xmu_1st_setup1 ; xmu_1st_setup2]);

figure
sgtitle('Mean of 18000 samples from posterior')
subplot(2,2,1)
showDistribution(xmu_0th_setup1,gridinfo); caxis(mucaxis); title('phi = [10 20 40], 0th order nn. Tikh')
subplot(2,2,2)
showDistribution(xmu_0th_setup2,gridinfo); caxis(mucaxis); title('phi = [10 20 70], 0th order nn. Tihk')
subplot(2,2,3)
showDistribution(xmu_1st_setup1,gridinfo); caxis(mucaxis); title('phi = [10 20 40], 1st order nn. Tikh')
subplot(2,2,4)
showDistribution(xmu_1st_setup2,gridinfo); caxis(mucaxis); title('phi = [10 20 70], 1st order nn. Tikh')


%% Show the four standard deviations with the same colorscale

stdcaxis(1) = min([xstd_0th_setup1 ; xstd_0th_setup2 ; xstd_1st_setup1 ; xstd_1st_setup2]);
stdcaxis(2) = max([xstd_0th_setup1 ; xstd_0th_setup2 ; xstd_1st_setup1 ; xstd_1st_setup2]);

figure
sgtitle('Standard deviation of 18000 samples from posterior')
subplot(2,2,1)
showDistribution(xstd_0th_setup1,gridinfo); caxis(stdcaxis); title('phi = [10 20 40], 0th order nn. Tikh')
subplot(2,2,2)
showDistribution(xstd_0th_setup2,gridinfo); caxis(stdcaxis); title('phi = [10 20 70], 0th order nn. Tihk')
subplot(2,2,3)
showDistribution(xstd_1st_setup1,gridinfo); caxis(stdcaxis); title('phi = [10 20 40], 1st order nn. Tikh')
subplot(2,2,4)
showDistribution(xstd_1st_setup2,gridinfo); caxis(stdcaxis); title('phi = [10 20 70], 1st order nn. Tikh')

















