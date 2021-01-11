clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_phi3')
    HPCDownload('biMax_sim_phi3','biMax_sim_phi3','s174483');
end

%%
%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('../data/biMax_sim_phi3')

files = {dir('../data/biMax_sim_phi3/*.mat').name}; 

%% Show regularisation parameters
clear, clc, close all
figure('units','normalized','outerposition',[0 0 1 1])
load('../data/biMax_sim_phi3/phi3_60_80_15.mat')
alphasim1 = alphasim1(1001:end);
qalpha1 = quantile(alphasim1,[0.025, 0.975]);
semilogx(alpha_relerr,r1st, 'LineWidth', 2, 'color', [0. 0.4470, 0.7410]); xlim([1e-2, 1e3]); hold on
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'o', 'LineWidth', 2, 'MarkerSize', 8, 'HandleVisibility','off', 'color', [0. 0.4470, 0.7410])
plot(alpha_relerr(idx1), min(r1st), '.', 'MarkerSize', 20, 'HandleVisibility', 'off', 'color', [0. 0.4470, 0.7410])
xline(qalpha1(1), '--', 'LineWidth', 2, 'HandleVisibility', 'off', 'color', [0. 0.4470, 0.7410])
xline(qalpha1(2), '--', 'LineWidth', 2, 'HandleVisibility', 'off', 'color', [0. 0.4470, 0.7410])

hold on; 
load('../data/biMax_sim_phi3/phi3_60_80_20.mat')
alphasim1 = alphasim1(1001:end);
qalpha1 = quantile(alphasim1,[0.025, 0.975]);
semilogx(alpha_relerr,r1st, 'LineWidth', 2, 'color', [0.8500, 0.3250, 0.0980]); xlim([1e-2, 1e3])
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'o', 'LineWidth', 2, 'MarkerSize', 8, 'HandleVisibility','off', 'color', [0.8500, 0.3250, 0.0980])
plot(alpha_relerr(idx1), min(r1st), '.', 'MarkerSize', 20, 'HandleVisibility', 'off', 'color', [0.8500, 0.3250, 0.0980])
xline(qalpha1(1), '--', 'LineWidth', 2, 'HandleVisibility', 'off', 'color', [0.8500, 0.3250, 0.0980])
xline(qalpha1(2), '--', 'LineWidth', 2, 'HandleVisibility', 'off', 'color', [0.8500, 0.3250, 0.0980])

load('../data/biMax_sim_phi3/phi3_60_80_25.mat')
alphasim1 = alphasim1(1001:end);
qalpha1 = quantile(alphasim1,[0.025, 0.975]);
semilogx(alpha_relerr,r1st, 'LineWidth', 2, 'color', [0.9290, 0.6940, 0.1250]); xlim([1e-2, 1e3])
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'o','LineWidth', 2,  'MarkerSize', 8, 'HandleVisibility','off', 'color', [0.9290, 0.6940, 0.1250])
plot(alpha_relerr(idx1), min(r1st), '.', 'MarkerSize', 20, 'HandleVisibility', 'off', 'color', [0.9290, 0.6940, 0.1250])
xline(qalpha1(1), '--','LineWidth', 2,  'HandleVisibility', 'off', 'color', [0.9290, 0.6940, 0.1250])
xline(qalpha1(2), '--','LineWidth', 2,  'HandleVisibility', 'off', 'color', [0.9290, 0.6940, 0.1250])

load('../data/biMax_sim_phi3/phi3_60_80_30.mat')
alphasim1 = alphasim1(1001:end);
qalpha1 = quantile(alphasim1,[0.025, 0.975]);
semilogx(alpha_relerr,r1st, 'LineWidth', 2, 'color', [0.4940, 0.1840, 0.5560]); xlim([1e-2, 1e3])
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'o', 'LineWidth', 2, 'MarkerSize', 8, 'HandleVisibility','off', 'color', 	[0.4940, 0.1840, 0.5560])
plot(alpha_relerr(idx1), min(r1st), '.', 'MarkerSize', 20, 'HandleVisibility', 'off', 'color', 	[0.4940, 0.1840, 0.5560])
xline(qalpha1(1), '--','LineWidth', 2,  'HandleVisibility', 'off', 'color', [0.4940, 0.1840, 0.5560])
xline(qalpha1(2), '--','LineWidth', 2,  'HandleVisibility', 'off', 'color', [0.4940, 0.1840, 0.5560])

ylim([0 0.35])
xlabel('\alpha','FontSize',15); ylabel('Relative error','FontSize',15)
legend('\phi_3 = 15','\phi_3 = 10','\phi_3 = 25','\phi_3 = 30','FontSize',15,'location','sw')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'biMax_phi3_analysis','-dpdf','-r0')


%% Check the convergence
clear, clc, close all
load('../data/biMax_sim_phi3/phi3_60_80_15.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

load('../data/biMax_sim_phi3/phi3_60_80_20.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

load('../data/biMax_sim_phi3/phi3_60_80_25.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

load('../data/biMax_sim_phi3/phi3_60_80_30.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

%% Check the convergence
clear, clc, close all
load('../data/biMax_sim_phi3/phi3_60_80_20.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

load('../data/biMax_sim_phi3/phi3_60_80_25.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

load('../data/biMax_sim_phi3/phi3_60_80_30.mat')
geweketest([alphasim1 deltasim1 lambdasim1], 1000)

%% Plot the mean reconstructions and credibility bounds
clear, clc, close all
load('../data/biMax_sim_phi3/phi3_60_80_15.mat', 'xNNHGS1');
xNNHGS1 = xNNHGS1(:,1001:end);
x15 = mean(xNNHGS1,2); xc15 = cbounds(xNNHGS1);

load('../data/biMax_sim_phi3/phi3_60_80_20.mat', 'xNNHGS1');
xNNHGS1 = xNNHGS1(:,1001:end);
x20 = mean(xNNHGS1,2); xc20 = cbounds(xNNHGS1);

load('../data/biMax_sim_phi3/phi3_60_80_25.mat', 'xNNHGS1');
xNNHGS1 = xNNHGS1(:,1001:end);
x25 = mean(xNNHGS1,2); xc25 = cbounds(xNNHGS1);

load('../data/biMax_sim_phi3/phi3_60_80_30.mat', 'xNNHGS1', 'ginfo');
xNNHGS1 = xNNHGS1(:,1001:end);
x30 = mean(xNNHGS1,2); xc30 = cbounds(xNNHGS1);
clear xNNHGS1

cax1 = [min([x15 ; x20 ; x25 ; x30]); max([x15 ; x20 ; x25 ; x30])];
cax2 = [min([xc15 ; xc20 ; xc25 ; xc30]); max([xc15 ; xc20 ; xc25 ; xc30])];

figure; showDistribution(x15,ginfo,cax1); print(gcf,'biMax_phi3_analysis_x15','-dpdf','-r0')
figure; showDistribution(xc15,ginfo,cax2); print(gcf,'biMax_phi3_analysis_xc15','-dpdf','-r0')

figure; showDistribution(x20,ginfo,cax1); print(gcf,'biMax_phi3_analysis_x20','-dpdf','-r0')
figure; showDistribution(xc20,ginfo,cax2); print(gcf,'biMax_phi3_analysis_xc20','-dpdf','-r0')

figure; showDistribution(x25,ginfo,cax1); print(gcf,'biMax_phi3_analysis_x25','-dpdf','-r0')
figure; showDistribution(xc25,ginfo,cax2); print(gcf,'biMax_phi3_analysis_xc25','-dpdf','-r0')

figure; showDistribution(x30,ginfo,cax1); print(gcf,'biMax_phi3_analysis_x30','-dpdf','-r0')
figure; showDistribution(xc30,ginfo,cax2); print(gcf,'biMax_phi3_analysis_xc30','-dpdf','-r0')


%% Make the uncertainty plot
clear, clc, close all
load('../data/biMax_sim_phi2/phi2_60_80_full.mat'); 
fprintf('%.3f & %.3f\n',min(r1st),relerr(xtrue,mean(xNNHGS1(:,501:end),2)))
load('../data/biMax_sim_phi3/phi3_60_80_15.mat'); 
fprintf('%.3f & %.3f\n',min(r1st),relerr(xtrue,mean(xNNHGS1(:,1001:end),2)))
xsim_60_80_15 = xNNHGS1(:,1001:end);
load('../data/biMax_sim_phi3/phi3_60_80_20.mat'); 
fprintf('%.3f & %.3f\n',min(r1st),relerr(xtrue,mean(xNNHGS1(:,1001:end),2)))
xsim_60_80_20 = xNNHGS1(:,1001:end);
load('../data/biMax_sim_phi3/phi3_60_80_25.mat'); 
fprintf('%.3f & %.3f\n',min(r1st),relerr(xtrue,mean(xNNHGS1(:,1001:end),2)))
xsim_60_80_25 = xNNHGS1(:,1001:end);
load('../data/biMax_sim_phi3/phi3_60_80_30.mat'); 
fprintf('%.3f & %.3f\n',min(r1st),relerr(xtrue,mean(xNNHGS1(:,1001:end),2)))
xsim_60_80_30 = xNNHGS1(:,1001:end);


xsim_60_80 = xNNHGS1(:,501:end);


%%
clear, clc, close all
load('../data/biMax_sim_phi2/phi2_60_80_full.mat','xNNHGS1','ginfo'); xbaseline = xNNHGS1(:,501:end);
load('../data/biMax_sim_phi3/phi3_60_80_15.mat', 'xNNHGS1'); x15 = xNNHGS1(:,1001:end);
load('../data/biMax_sim_phi3/phi3_60_80_20.mat', 'xNNHGS1'); x20 = xNNHGS1(:,1001:end);
load('../data/biMax_sim_phi3/phi3_60_80_25.mat', 'xNNHGS1'); x25 = xNNHGS1(:,1001:end);
load('../data/biMax_sim_phi3/phi3_60_80_30.mat', 'xNNHGS1'); x30 = xNNHGS1(:,1001:end);

figure
UQmap(x15,xbaseline,ginfo); caxis([-2e5, 2e5])

figure
UQmap(x30,xbaseline,ginfo); caxis([-2e5, 2e5])

%% Construct the analysis tables
phi_analysis_table(files)

%% Show the different plots
phi_analysis_plot(files)

