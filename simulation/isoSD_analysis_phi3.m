clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

files = {dir('./data/isoSD_sim_phi3/*.mat').name}; 

[T0, T1] = phi_analysis_table(files); close all

%phi_analysis_plot(files)

T0new = T0;
T0new(:,2:3) = [];
T0new(:,end-2:end) = [];
T0new

T1new = T1;
T1new(:,2:3) = []; T1new(:,end-2:end) = [];
T1new

table2latex(T0new,'./data/phi3T0_isoSD.tex')
table2latex(T1new,'./data/phi3T1_isoSD.tex')

%% 
clear, clc, close all
load('phi3_60_80_75.mat')
    
figure
semilogx(alpha_relerr, r1st, 'b'); xlim([1e0, 1e10])
hold on
plot(optalpha_1st, minr1, 'b.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'bo', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim1),relerr(xtrue, xsamplealpha1), ...
         qalpha1(1)-mean(alphasim1), qalpha1(2)-mean(alphasim1), ...
         'horizontal','LineWidth',1,'color','b', 'HandleVisibility', 'off')

%%
clear, clc, close all
load('phi3_60_80_15.mat')
xsim_60_80_15 = xNNHGS1;

load('phi3_60_80_30.mat')
xsim_60_80_30 = xNNHGS1;

figure
UQmap(xsim_60_80_15, xsim_60_80_30, ginfo, 'cbound', false)

%% Show the different slowing down distributions
clear, clc, close all

load('phi3_60_80_15.mat')
xsim_60_80_15 = xNNHGS1; 

load('phi3_60_80_20.mat')
xsim_60_80_20 = xNNHGS1; 

load('phi3_60_80_25.mat')
xsim_60_80_25 = xNNHGS1; 

load('phi3_60_80_30.mat')
xsim_60_80_30 = xNNHGS1; 

cax = [min([xsim_60_80_15(:,2);xsim_60_80_20(:,2);xsim_60_80_25(:,2);xsim_60_80_30(:,2)]) ...
       max([xsim_60_80_15(:,2);xsim_60_80_20(:,2);xsim_60_80_25(:,2);xsim_60_80_30(:,2)])];
       
figure
showDistribution(xsim_60_80_15(:,2),ginfo, cax)

figure
showDistribution(xsim_60_80_20(:,2),ginfo, cax)

figure
showDistribution(xsim_60_80_25(:,2),ginfo, cax)

figure
showDistribution(xsim_60_80_30(:,2),ginfo, cax)

