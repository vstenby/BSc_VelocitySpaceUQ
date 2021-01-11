%% Analysis of the data sim folder.

clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('../data/data_sim')
    HPCDownload('data_sim','./data/data_sim','s174483');
end

% Set the timestamp for analysis
ts = 4;
[sim_rhs_01, sim_rhs_02, sim_rhs_03, sim_rhs_04] = load_data(ts, '../data/data_sim');

%%

figure
semilogx(sim_rhs_04.alpha, sim_rhs_04.relerr, 'k-')
hold on
qalpha(1) = quantile(alphasim,0.025);
qalpha(2) = mean(alphasim);
qalpha(3) = quantile(alphasim,0.975);
xline(qalpha(1), 'b-')
xline(qalpha(2), 'b-')
xline(qalpha(3), 'b-')

qalpha(1) = quantile(sim_rhs_04.alphasim,0.025);
qalpha(2) = mean(sim_rhs_04.alphasim);
qalpha(3) = quantile(sim_rhs_04.alphasim,0.975);
xline(qalpha(1), 'r-')
xline(qalpha(2), 'r-')
xline(qalpha(3), 'r-')


figure
subplot(1,2,1)
showDistribution(xsim(:,1))
subplot(1,2,2)
showDistribution(sim_rhs_03.xsim(:,1))
%% 

figure
showDistribution(sim_rhs_01.xtrue); 
%saveas(gcf, 'ts_06_xtrue.pdf')

%%
figure
plot(sim_rhs_01.b);
hold on
plot(sim_rhs_02.b); plot(sim_rhs_03.b); plot(sim_rhs_04.b);
title('Right hand sides')
legend('Inverse crime', 'No inverse crime', 'FIDAsim', 'Real data', 'location', 'northwest','FontSize',15)
%saveas(gcf, 'ts_06_rhs.pdf')

%% Various reconstructions
close all
[~,min_rhs1] = min(sim_rhs_01.relerr);
[~,min_rhs2] = min(sim_rhs_02.relerr);
[~,min_rhs3] = min(sim_rhs_03.relerr);
[~,min_rhs4] = min(sim_rhs_04.relerr);

switch ts
    case 6
        idx_invcrime = [1, 20, 40, min_rhs1, 100, 150];
        idx_nocrime  = [1, 20, 35, min_rhs2, 130, 170];
        idx_fidasim  = [1, 40, 50, min_rhs3, 150, 190];
        idx_realdata = [1, 20, 35, 75, min_rhs4, 140];
    otherwise     
end

figure('units','normalized','outerposition',[0 0 1 1])

recon_plot(sim_rhs_01, idx_invcrime); sgtitle('Inverse crime','FontSize',20)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_recon_invcrime','-dpdf','-r0')

figure('units','normalized','outerposition',[0 0 1 1])
recon_plot(sim_rhs_02, idx_nocrime); sgtitle('No inverse crime','FontSize',20)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_recon_nocrime','-dpdf','-r0')

figure('units','normalized','outerposition',[0 0 1 1])
recon_plot(sim_rhs_03, idx_fidasim); sgtitle('FIDAsim','FontSize',20)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_recon_fidasim','-dpdf','-r0')

figure('units','normalized','outerposition',[0 0 1 1])
recon_plot(sim_rhs_04, idx_realdata); sgtitle('Real data','FontSize',20)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_recon_realdata','-dpdf','-r0')

%% Chain plots

chain_analysis(sim_rhs_01.deltasim, sim_rhs_01.lambdasim)
chain_analysis(sim_rhs_02.deltasim, sim_rhs_02.lambdasim)
chain_analysis(sim_rhs_03.deltasim, sim_rhs_03.lambdasim)
chain_analysis(sim_rhs_04.deltasim, sim_rhs_04.lambdasim)

%%
figure
chain_plot(sim_rhs_01, 'Inverse crime')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_chainplot_invcrime','-dpdf','-r0')

figure
chain_plot(sim_rhs_02, 'No inverse crime')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_chainplot_nocrime','-dpdf','-r0')

figure
chain_plot(sim_rhs_03, 'FIDAsim')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_chainplot_fidasim','-dpdf','-r0')

figure
chain_plot(sim_rhs_04, 'Real data')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_04_chainplot_realdata','-dpdf','-r0')

%%

figure('units','normalized','outerposition',[0 0 1 1])
NNHGS_plot(sim_rhs_01); 
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_06_UQ_invcrime','-dpdf','-r0')

figure('units','normalized','outerposition',[0 0 1 1])
NNHGS_plot(sim_rhs_02); 
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_06_UQ_nocrime','-dpdf','-r0')

figure('units','normalized','outerposition',[0 0 1 1])
NNHGS_plot(sim_rhs_03); 
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_06_UQ_fidasim','-dpdf','-r0')

figure('units','normalized','outerposition',[0 0 1 1])
NNHGS_plot(sim_rhs_04); 
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'ts_06_UQ_realdata','-dpdf','-r0')

%% Auxil functions

function [sim_rhs_01, sim_rhs_02, sim_rhs_03, sim_rhs_04] = load_data(ts, folder)
    %Load the simulation data.
    for rhs = 1:4
        load(sprintf('./%s/ts_%02d_rhs_%02d.mat', folder, ts, rhs))

        %Load the true solution
        eval(sprintf('sim_rhs_%02d.xtrue = xtrue(:);',rhs))

        %Load the Tikhonov results
        eval(sprintf('sim_rhs_%02d.alpha = alpha;',rhs))
        eval(sprintf('sim_rhs_%02d.xalpha = xalpha;',rhs))
        eval(sprintf('sim_rhs_%02d.relerr = relerr;',rhs))

        %Load the NNHGS-results.
        eval(sprintf('sim_rhs_%02d.xsim = xsim;',rhs))
        eval(sprintf('sim_rhs_%02d.alphasim = alphasim;',rhs))
        eval(sprintf('sim_rhs_%02d.deltasim = deltasim;',rhs))
        eval(sprintf('sim_rhs_%02d.lambdasim = lambdasim;',rhs))
        eval(sprintf('sim_rhs_%02d.info = info;',rhs))

        clearvars -except sim_rhs* ts rhs folder
    end

    %Note that in the case of rhs_04, A and b is scaled by the error.
    load(sprintf('../data/TCV/ts%d.mat',ts))
    sim_rhs_01.b = double(S_TR_invcrime_n);
    sim_rhs_02.b = double(S_TR_n);
    sim_rhs_03.b = double(fida_total_n);
    sim_rhs_04.b = double(S_data);

    clearvars -except sim_rhs* 
end

function showDistribution(x)
imagesc(reshape(x,20,20)); axis image; axis xy; colorbar()
end

function recon_plot(sim, idx)
colors = [128 0 0; ...
          128 128 0; ...
          0 128 128; ...
          0 0 128; ...
          230 25 75; ...
          245 130 48]./255;
cax = [min([sim.xtrue ; sim.xalpha(:,idx(1)) ; sim.xalpha(:,idx(2)) ; sim.xalpha(:,idx(3)) ; ...
            sim.xalpha(:,idx(4)) ; sim.xalpha(:,idx(5)); sim.xalpha(:,idx(6))]); ...
       max([sim.xtrue ; sim.xalpha(:,idx(1)) ; sim.xalpha(:,idx(2)) ; sim.xalpha(:,idx(3)) ; ...
            sim.xalpha(:,idx(4)) ; sim.xalpha(:,idx(5)); sim.xalpha(:,idx(6))])]; 
%Show the different reconstructions
subplot(3,4,[5,6,9,10])
semilogx(sim.alpha, sim.relerr, 'k-'); title('\alpha vs relative error')
hold on
scatter(sim.alpha(idx), sim.relerr(idx), 50, colors, 'filled')
subplot(3,4,[1,2])
showDistribution(sim.xtrue); title('True solution');
subplot(3,4,3)
showDistribution(sim.xalpha(:,idx(1))); title(sprintf('alpha = %.3e, r = %.2f',sim.alpha(idx(1)), sim.relerr(idx(1))),'color',colors(1,:))
subplot(3,4,4)
showDistribution(sim.xalpha(:,idx(2))); title(sprintf('alpha = %.3e, r = %.2f',sim.alpha(idx(2)), sim.relerr(idx(2))),'color',colors(2,:))
subplot(3,4,7)
showDistribution(sim.xalpha(:,idx(3))); title(sprintf('alpha = %.3e, r = %.2f',sim.alpha(idx(3)), sim.relerr(idx(3))),'color',colors(3,:))
subplot(3,4,8)
showDistribution(sim.xalpha(:,idx(4))); title(sprintf('alpha = %.3e, r = %.2f',sim.alpha(idx(4)), sim.relerr(idx(4))),'color',colors(4,:))
subplot(3,4,11)
showDistribution(sim.xalpha(:,idx(5))); title(sprintf('alpha = %.3e, r = %.2f',sim.alpha(idx(5)), sim.relerr(idx(5))),'color',colors(5,:))
subplot(3,4,12)
showDistribution(sim.xalpha(:,idx(6))); title(sprintf('alpha = %.3e, r = %.2f',sim.alpha(idx(6)), sim.relerr(idx(6))),'color',colors(6,:))
end

function chain_plot(sim,s)
subplot(2,3,1)
plot(sim.alphasim)
xlim([1 length(sim.alphasim)]); title('\alpha')
xtickangle(45)

subplot(2,3,2)
plot(sim.deltasim)
xlim([1 length(sim.deltasim)]); title('\delta')
xtickangle(45)

subplot(2,3,3)
plot(sim.lambdasim)
xlim([1 length(sim.lambdasim)]); title('\lambda')
xtickangle(45)

subplot(2,3,4)
hist(sim.alphasim); title('\alpha')

subplot(2,3,5)
hist(sim.deltasim); title('\delta')

subplot(2,3,6)
hist(sim.lambdasim); title('\lambda')

sgtitle(sprintf('%s | Geweke p = %f',s, geweketest([sim.alphasim sim.deltasim sim.lambdasim])))
end

function NNHGS_plot(sim)

%Find the minimum reconstruction
[~,idx] = min(sim.relerr);
qalpha    = quantile(sim.alphasim, [0.025, 0.975]);
meanalpha = mean(sim.alphasim);

%Show the different reconstructions
subplot(2,3,[1,2])
semilogx(sim.alpha, sim.relerr, 'k-'); title('\alpha vs relative error','FontSize',30)
hold on
plot(sim.alpha(idx),sim.relerr(idx),'k.','MarkerSize',15)
xline(qalpha(1), 'k--')
xline(qalpha(2), 'k--')
xline(meanalpha, 'k-')

%cax = [min([sim.xtrue ; sim.xalpha(:,idx) ; sim.xsim(:,1)])
subplot(2,3,3)
showDistribution(sim.xtrue); title('True solution', 'FontSize', 30)

subplot(2,3,4)
showDistribution(sim.xalpha(:,idx)); title(sprintf('Minimum r solution, r = %.2f',sim.relerr(idx)), 'FontSize', 30)

%Posterior mean
subplot(2,3,5)
showDistribution(sim.xsim(:,1)); title(sprintf('Sample mean, r = %.2f', relerr(sim.xtrue, sim.xsim(:,1))), 'FontSize', 30)

%Standard deviation
subplot(2,3,6)
showDistribution(sim.xsim(:,2)); title('Standard deviation', 'FontSize', 30)

end


