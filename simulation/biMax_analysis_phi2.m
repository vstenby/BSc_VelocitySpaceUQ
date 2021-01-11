%% Analysis of the phi2 experiment.
clear, clc, close all

addpath('./data/biMax_sim_phi2')

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Grab the relevant simulations
rownames = combvec(10:10:80,10:10:80)';
d = unique(sort(rownames,2), 'rows');
d(d(:,1)==d(:,2),:) = [];

files = cell(size(d,1),1);
for i=1:size(d,1)
     files{i} = sprintf('./data/biMax_sim_phi2/phi2_%02d_%02d.mat',d(i,1),d(i,2));
end
%% Construct the analysis tables
[T0, T1] = phi_analysis_table(files);
close all

T0new = T0;
T0new(:,1:3) = [];

%table2latex(T0new, './data/phi2T0.tex')

%table2latex(T0new,'phi2T0.tex')

%% Convergence plots for one of the simulations
clear, clc, close all

load('./data/biMax_sim_phi2/phi2_60_80.mat')

chain_analysis(lambdasim0, deltasim0)

%%

T1new = T1;
T1new(:,1:3) = [];

table2latex(T1new, './data/phi2T1.tex')

%% Count cases
[status, ~, group] = unique(T1new(:,5).Variables);

sprintf('%s : %d', status{1}, sum(group==1))
sprintf('%s : %d', status{2}, sum(group==2))
sprintf('%s : %d', status{3}, sum(group==3))

%% 0th Tikhonov Investigation
figure
load('phi2_10_40.mat','alpha_relerr','r0th', 'optalpha_0th', 'minr0', 'alphasim0', 'xtrue', 'xsamplealpha0', 'qalpha0')
semilogx(alpha_relerr,r0th,'r-'); xlim([1e-20, 1e-9])
hold on
plot(optalpha_0th, minr0, 'r.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim0), relerr(xtrue, xsamplealpha0), 'ro', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim0),relerr(xtrue, xsamplealpha0), ...
         qalpha0(1)-mean(alphasim0), qalpha0(2)-mean(alphasim0), ...
         'horizontal','LineWidth',1,'color','r', 'HandleVisibility', 'off')

load('phi2_20_50.mat','r0th', 'optalpha_0th', 'minr0', 'alphasim0', 'xtrue', 'xsamplealpha0', 'qalpha0')
semilogx(alpha_relerr,r0th, 'b-'); 
hold on
plot(optalpha_0th, minr0, 'b.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim0), relerr(xtrue, xsamplealpha0), 'bo', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim0),relerr(xtrue, xsamplealpha0), ...
         qalpha0(1)-mean(alphasim0), qalpha0(2)-mean(alphasim0), ...
         'horizontal','LineWidth',1,'color','b', 'HandleVisibility', 'off')

load('phi2_60_80.mat','r0th', 'optalpha_0th', 'minr0', 'alphasim0', 'xtrue', 'xsamplealpha0', 'qalpha0')
semilogx(alpha_relerr,r0th,'k-')
hold on
plot(optalpha_0th, minr0, 'k.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim0), relerr(xtrue, xsamplealpha0), 'ko', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim0),relerr(xtrue, xsamplealpha0), ...
         qalpha0(1)-mean(alphasim0), qalpha0(2)-mean(alphasim0), ...
         'horizontal','LineWidth',1,'color','k', 'HandleVisibility', 'off')


legend('10 40', '20 50', '60 80', 'location', 'nw','FontSize',15)
xlabel('\alpha','FontSize',20); ylabel('Relative error','FontSize',20)

%% 1st Tikhonov Investigation
figure
load('phi2_60_80.mat','alpha_relerr','r1st', 'optalpha_1st', 'minr1', 'alphasim1', 'xtrue', 'xsamplealpha1', 'qalpha1')
semilogx(alpha_relerr,r1st,'r-'); xlim([1e-4, 1e5])
hold on
plot(optalpha_1st, minr1, 'r.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'ro', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim1),relerr(xtrue, xsamplealpha1), ...
         qalpha1(1)-mean(alphasim1), qalpha1(2)-mean(alphasim1), ...
         'horizontal','LineWidth',1,'color','r', 'HandleVisibility', 'off')

load('phi2_50_80.mat','r1st', 'optalpha_1st', 'minr1', 'alphasim1', 'xtrue', 'xsamplealpha1', 'qalpha1')
semilogx(alpha_relerr,r1st, 'b-'); 
hold on
plot(optalpha_1st, minr1, 'b.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'bo', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim1),relerr(xtrue, xsamplealpha1), ...
         qalpha1(1)-mean(alphasim1), qalpha1(2)-mean(alphasim1), ...
         'horizontal','LineWidth',1,'color','b', 'HandleVisibility', 'off')

load('phi2_20_40.mat','r1st', 'optalpha_1st', 'minr1', 'alphasim1', 'xtrue', 'xsamplealpha1', 'qalpha1')
semilogx(alpha_relerr,r1st,'k-')
hold on
plot(optalpha_1st, minr1, 'k.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim1), relerr(xtrue, xsamplealpha1), 'ko', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim1),relerr(xtrue, xsamplealpha1), ...
         qalpha1(1)-mean(alphasim1), qalpha1(2)-mean(alphasim1), ...
         'horizontal','LineWidth',1,'color','k', 'HandleVisibility', 'off')


legend('60 80 (Too much)', '50 80 (Good)', '20 40 (Too little)', 'location', 'nw','FontSize',15)
xlabel('\alpha','FontSize',20); ylabel('Relative error','FontSize',20)
%% Show x 
clear, clc, close all

load('phi2_10_50.mat', 'xNNHGS1', 'ginfo', 'xtrue'); x_10_50 = xNNHGS1(:,1);
load('phi2_50_80.mat', 'xNNHGS1', 'ginfo'); x_50_80 = xNNHGS1(:,1);
load('phi2_20_50.mat', 'xNNHGS1', 'ginfo'); x_20_50 = xNNHGS1(:,1);
load('phi2_30_80.mat', 'xNNHGS1', 'ginfo'); x_30_80 = xNNHGS1(:,1);
load('phi2_10_30.mat', 'xNNHGS1', 'ginfo'); x_10_30 = xNNHGS1(:,1);
load('phi2_60_80.mat', 'xNNHGS1', 'ginfo'); x_60_80 = xNNHGS1(:,1);

caxis = [min([x_10_50;x_50_80;x_20_50;x_30_80;x_10_30;x_60_80]), ...
         max([x_10_50;x_50_80;x_20_50;x_30_80;x_10_30;x_60_80])];
figure
showDistribution(x_10_50, ginfo, caxis);

figure
showDistribution(x_50_80, ginfo, caxis);

figure
showDistribution(x_20_50, ginfo, caxis);

figure
showDistribution(x_30_80, ginfo, caxis);

figure
showDistribution(x_10_30, ginfo, caxis);

figure
showDistribution(x_60_80, ginfo, caxis);

figure
showDistribution(xtrue, ginfo, caxis);
%% Further analysis of phi2_10_30_full
clear, clc, close all
HPCDownload('./biMax_sim_phi2/phi2_10_50_full.mat','./data/biMax_sim_phi2/phi2_10_50_full.mat','s174483')
HPCDownload('./biMax_sim_phi2/phi2_50_80_full.mat','./data/biMax_sim_phi2/phi2_50_80_full.mat','s174483')
HPCDownload('./biMax_sim_phi2/phi2_20_50_full.mat','./data/biMax_sim_phi2/phi2_20_50_full.mat','s174483')
HPCDownload('./biMax_sim_phi2/phi2_30_80_full.mat','./data/biMax_sim_phi2/phi2_30_80_full.mat','s174483')
HPCDownload('./biMax_sim_phi2/phi2_10_30_full.mat','./data/biMax_sim_phi2/phi2_10_30_full.mat','s174483')
HPCDownload('./biMax_sim_phi2/phi2_60_80_full.mat','./data/biMax_sim_phi2/phi2_60_80_full.mat','s174483')

%%
clear, clc, close all
phis = {'phi2_10_50_full.mat','phi2_50_80_full.mat','phi2_20_50_full.mat',...
        'phi2_30_80_full.mat','phi2_10_30_full.mat','phi2_60_80_full.mat'};
    
qwidth = zeros(450,6);
for i=1:length(phis)
    load(phis{i}, 'xNNHGS1', 'xtrue', 'ginfo')
    xsim1 = xNNHGS1(:,501:end);

    qlower = quantile(xsim1,0.05/2,2);
    qupper = quantile(xsim1,1-(0.05/2),2);
    qwidth(:, i) = qupper - qlower;
end

caxis = [min(qwidth(:)), max(qwidth(:))];

for i=1:length(phis)
   figure(i)
   showDistribution(qwidth(:,i),ginfo,caxis)
end

%%

clear, clc, close all
load('phi2_10_50_full.mat', 'xNNHGS1')
xsim_10_50_1st = xNNHGS1(:,501:end);

load('phi2_60_80_full.mat', 'xNNHGS1', 'ginfo')
xsim_60_80_1st = xNNHGS1(:,501:end);

UQmap(xsim_10_50_1st, xsim_60_80_1st, ginfo)
%%
clear, clc, close all
load('phi2_10_50_full.mat', 'xNNHGS1')
xsim1 = xNNHGS1(:,501:end);

qlower = quantile(xsim1,0.05/2,2);
qupper = quantile(xsim1,1-(0.05/2),2);
qwidth_10_50 = qupper - qlower;

load('phi2_60_80_full.mat', 'xNNHGS1', 'ginfo')
xsim1 = xNNHGS1(:,501:end);

qlower = quantile(xsim1,0.05/2,2);
qupper = quantile(xsim1,1-(0.05/2),2);
qwidth_60_80 = qupper - qlower;

dif_uncertainty = qwidth_60_80-qwidth_10_50;
%if dif_uncertainty > 0, then 10 50 is more certain. 
%if dif_uncertainty < 0, then 60 80 is more certain.



showDistribution(dif_uncertainty,ginfo); 
c = colorbar;
c.Ticks = [min(dif_uncertainty), 0, max(dif_uncertainty)];
%c.TickLabels = {'60 80', 'Same', '10 50'};
%colormap('hot')
%ylabel(c, 'Lowest uncertainty in x')

%% Which simulation has the highest uncertainty?
clear, clc, close all
load('phi2_10_50_full.mat', 'xNNHGS1', 'ginfo'); xNNHGS1_10_50 = xNNHGS1(:,501:end);
load('phi2_50_80_full.mat', 'xNNHGS1'); xNNHGS1_50_80 = xNNHGS1(:,501:end);
load('phi2_20_50_full.mat', 'xNNHGS1'); xNNHGS1_20_50 = xNNHGS1(:,501:end);
load('phi2_30_80_full.mat', 'xNNHGS1'); xNNHGS1_30_80 = xNNHGS1(:,501:end);
load('phi2_10_30_full.mat', 'xNNHGS1'); xNNHGS1_10_30 = xNNHGS1(:,501:end);
load('phi2_60_80_full.mat', 'xNNHGS1'); xNNHGS1_60_80 = xNNHGS1(:,501:end);


[~, I] = min([cbounds(xNNHGS1_10_50) cbounds(xNNHGS1_50_80) ...
              cbounds(xNNHGS1_20_50) cbounds(xNNHGS1_30_80) ...
              cbounds(xNNHGS1_10_30) cbounds(xNNHGS1_60_80)],[],2);

showDistribution(I,ginfo)
c = colorbar;
c.TickLabels = {'10 50', '50 80', '20 50', '30 80', '10 30', '60 80'};

%% Do the same analysis for slowing down distribution
clear, clc, close all
HPCDownload('./isoSD_sim_phi2/phi2_10_50.mat', './data/isoSD_sim_phi2/phi2_10_50.mat', 's174483')
HPCDownload('./isoSD_sim_phi2/phi2_60_80.mat', './data/isoSD_sim_phi2/phi2_60_80.mat', 's174483')

%%
clear, clc, close all
load('./data/isoSD_sim_phi2/phi2_10_50.mat')

cbounds_10_50 = cbounds(xNNHGS0(:,501:end));

load('./data/isoSD_sim_phi2/phi2_60_80.mat')

cbounds_60_80 = cbounds(xNNHGS0(:,501:end));

dif_uncertainty = cbounds_60_80-cbounds_10_50;
showDistribution(dif_uncertainty,ginfo); 
c = colorbar;
c.Ticks = [min(dif_uncertainty), 0, max(dif_uncertainty)];
c.TickLabels = {'60 80', 'Same', '10 50'};
colormap('hot')
ylabel(c, 'Model with lowest uncertainty')



%%

function c = cbounds(xsamp)
% npixel x nsamps matrix
qlower = quantile(xsamp,0.05/2,2);
qupper = quantile(xsamp,1-(0.05/2),2);
c = qupper - qlower;
end







