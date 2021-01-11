clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))
addpath(genpath('../data'))

gewekeps0 = zeros(6,1);
alphaopt0 = zeros(6,1);
alphamean0 = zeros(6,1);
qalpha0_1 = zeros(6,1);
qalpha0_2 = zeros(6,1);
status0 = cell(6,1);

for i=1:6
    load(sprintf('DentalCT_setup%02d.mat',i))

    gewekeps0(i) = geweketest([alphasim0(501:end), deltasim0(501:end), lambdasim0(501:end)]);
    alphaopt0(i) = optalpha_0th;
    alphamean0(i) = mean(alphasim0(501:end));
    qalpha0_1(i) = qalpha0(1);
    qalpha0_2(i) = qalpha0(2);

    if alphaopt0(i) < qalpha0_1(i)
       %Overregularized
       status0{i} = 'Too much';
    elseif alphaopt0(i) > qalpha0_2(i)
       status0{i} = 'Too little';
    else
       status0{i} = 'Good';
    end
end

setups = {'Narrow left';'Wide left';'Full';'Front';'Wide right';'Narrow right'};

T0 = table(setups, gewekeps0, alphaopt0, alphamean0, qalpha0_1, qalpha0_2, status0);

table2latex(T0,'./data/dentalCT_T0.tex')
%%

gewekeps1 = zeros(6,1);
alphaopt1 = zeros(6,1);
alphamean1 = zeros(6,1);
qalpha1_1 = zeros(6,1);
qalpha1_2 = zeros(6,1);
status1 = cell(6,1);

for i=1:6
    load(sprintf('DentalCT_setup%02d.mat',i))

    gewekeps1(i) = geweketest([alphasim1(501:end), deltasim1(501:end), lambdasim1(501:end)]);
    alphaopt1(i) = optalpha_1st;
    alphamean1(i) = mean(alphasim1(501:end));
    qalpha1_1(i) = qalpha1(1);
    qalpha1_2(i) = qalpha1(2);

    if alphaopt1(i) < qalpha1_1(i)
       %Overregularized
       status1{i} = 'Too much';
    elseif alphaopt1(i) > qalpha1_2(i)
       status1{i} = 'Too little';
    else
       status0{i} = 'Good';
    end
end

setups = {'Narrow left';'Wide left';'Full';'Front';'Wide right';'Narrow right'};

T1 = table(setups, gewekeps1, alphaopt1, alphamean1, qalpha1_1, qalpha1_2, status1);

table2latex(T1,'./data/dentalCT_T1.tex')

%% Relative errors
clear, clc, close all

minr0s = zeros(6,1); r_x_alphamean0 = zeros(6,1); r_xmean0 = zeros(6,1);
minr1s = zeros(6,1); r_x_alphamean1 = zeros(6,1); r_xmean1 = zeros(6,1);

for i=1:6
    load(sprintf('DentalCT_setup%02d.mat',i))
    minr0s(i) = minr0; minr1s(i) = minr1;
    r_x_alphamean0(i) = relerr(x, xsamplealpha0);
    r_x_alphamean1(i) = relerr(x, xsamplealpha1);
    r_xmean0(i) = relerr(x, mean(xNNHGS0(:,501:end),2));
    r_xmean1(i) = relerr(x, mean(xNNHGS1(:,501:end),2));
end

T = table(minr0s, r_x_alphamean0, r_xmean0, minr1s, r_x_alphamean1, r_xmean1);
T.Properties.RowNames = {'Narrow left', 'Wide left', 'Full', 'Front', 'Wide right', 'Narrow right'};
table2latex(T, './data/DentalCT_r.tex')

%% 0th order investigation
%Plot the chains
nburnin = 500;
chain_analysis(deltasim0(501:end),lambdasim0(501:end))

figure
semilogx(alphavec, r0, 'r'); 
hold on
plot(optalpha_0th, minr0, 'r.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim0), relerr(x, xsamplealpha0), 'ro', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim0),relerr(x, xsamplealpha0), ...
         qalpha0(1)-mean(alphasim0), qalpha0(2)-mean(alphasim0), ...
         'horizontal','LineWidth',1,'color','r', 'HandleVisibility', 'off')

%%
figure
imagesc(reshape(xopt0,30,30)); axis image; title(sprintf('relerr = %.2f',minr0))

figure
imagesc(reshape(xsamplealpha0,30,30)); axis image; title(sprintf('relerr = %.2f',relerr(x,xsamplealpha0)))

figure
imagesc(reshape(mean(xNNHGS0(:,501:end),2),30,30)); axis image; title(sprintf('relerr = %.2f',relerr(x,mean(xNNHGS0(:,501:end),2))))
     
qalpha = quantile(xNNHGS0(:,501:end),[0.025, 0.975],2); qwidth0 = qalpha(:,2)-qalpha(:,1);
figure
imagesc(reshape(qwidth0,30,30)); axis image; colorbar()

%% Save all of the plots!
clear, clc, close all
for i=1:6
    load(sprintf('DentalCT_setup%02d.mat',i))
    xmean0 = mean(xNNHGS0(:,501:end),2);
    xmean1 = mean(xNNHGS1(:,501:end),2);
    xopt0 = xopt0;
    xopt1 = xopt1;
    c0 = cbounds(xNNHGS0(:,501:end));
    c1 = cbounds(xNNHGS1(:,501:end));
    
    caxis1 = [min([xmean0 ; xmean1 ; xopt0 ; xopt1]), max([xmean0 ; xmean1 ; xopt0 ; xopt1])];
    caxis2 = [min([c0 ; c1]), max([c0 ; c1])];
    
    figure
    imagesc(reshape(xopt0,30,30)); axis image; caxis(caxis1); colorbar()
    saveas(gcf,sprintf('./data/dentalctplots/xopt0_setup%02d.pdf',i))
  
    figure
    imagesc(reshape(xmean0,30,30)); axis image; caxis(caxis1); colorbar()
    saveas(gcf,sprintf('./data/dentalctplots/xmean0_setup%02d.pdf',i))
    
    figure
    imagesc(reshape(c0,30,30)); axis image; caxis(caxis2); colorbar()
    saveas(gcf,sprintf('./data/dentalctplots/xcbounds0_setup%02d.pdf',i))
    
    figure
    imagesc(reshape(xopt1,30,30)); axis image; caxis(caxis1); colorbar()
    saveas(gcf,sprintf('./data/dentalctplots/xopt1_setup%02d.pdf',i))
    
    figure
    imagesc(reshape(xmean1,30,30)); axis image; caxis(caxis1); colorbar()
    saveas(gcf,sprintf('./data/dentalctplots/xmean1_setup%02d.pdf',i))
    
    figure
    imagesc(reshape(c1,30,30)); axis image; caxis(caxis2); colorbar()
    saveas(gcf,sprintf('./data/dentalctplots/xcbounds1_setup%02d.pdf',i))
    clc
end

%% Compare the uncertainties
clear, clc, close all
load('./data/DentalCT_sim/DentalCT_setup01.mat')
c01_0th = cbounds(xNNHGS0(:,501:end));
c01_1st = cbounds(xNNHGS1(:,501:end));

load('./data/DentalCT_sim/DentalCT_setup03.mat')
c03_0th = cbounds(xNNHGS0(:,501:end));
c03_1st = cbounds(xNNHGS1(:,501:end));


figure
dif = c01_1st-c03_1st; 
imagesc(reshape(dif,30,30)); axis image; colorbar()

hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')

caxis([-max(abs(dif)), max(abs(dif))])
colormap(brewermap([],'*RdBu'))

%%
clear, clc, close all
load('../data/DentalCT_sim/DentalCT_setup03.mat','xNNHGS1');
xbaseline = xNNHGS1(:,501:end);

load('../data/DentalCT_sim/DentalCT_setup01.mat','xNNHGS1');
x01_1st = xNNHGS1(:,501:end);
clear xxNHGS1

load('../data/DentalCT_sim/DentalCT_setup02.mat','xNNHGS1');
x02_1st = xNNHGS1(:,501:end);
clear xxNHGS1

load('../data/DentalCT_sim/DentalCT_setup04.mat','xNNHGS1');
x04_1st = xNNHGS1(:,501:end);
clear xxNHGS1

load('../data/DentalCT_sim/DentalCT_setup05.mat','xNNHGS1');
x05_1st = xNNHGS1(:,501:end);
clear xxNHGS1

load('../data/DentalCT_sim/DentalCT_setup06.mat','xNNHGS1');
x06_1st = xNNHGS1(:,501:end);
clear xxNHGS1

figure
UQmap(x01_1st,xbaseline,[30,30],'display','imagesc')
hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')


figure
UQmap(x02_1st,xbaseline,[30,30],'display','imagesc')
hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')

figure
UQmap(x04_1st,xbaseline,[30,30],'display','imagesc')
hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')


figure
UQmap(x05_1st,xbaseline,[30,30],'display','imagesc')
hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')


figure
UQmap(x06_1st,xbaseline,[30,30],'display','imagesc')
hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')


%%
imagesc(reshape(cbounds(x01_1st) - cbounds(x03_1st),30,30)); axis image; colorbar()



%%
clear, clc, close all
load('./data/DentalCT_sim/DentalCT_setup03.mat','xNNHGS0','xNNHGS1');




%%

clear, clc, close all
load('../data/DentalCT_sim/DentalCT_setup02.mat','xNNHGS0');
x01_1st = xNNHGS0(:,501:end);
load('../data/DentalCT_sim/DentalCT_setup03.mat','xNNHGS0');
x03_1st = xNNHGS0(:,501:end);
clear xxNHGS1

UQmap(x01_1st,x03_1st,[30,30])
hold on
[~,~,mask] = dental_phantom(30);
visboundaries(mask,'color','k')


%% 1st order investigation
figure
semilogx(alphavec, r1, 'b'); 
hold on
plot(optalpha_1st, minr1, 'b.', 'MarkerSize', 15, 'HandleVisibility','off')
plot(mean(alphasim1), relerr(x, xsamplealpha1), 'bo', 'MarkerSize', 8, 'HandleVisibility','off')
errorbar(mean(alphasim1),relerr(x, xsamplealpha1), ...
         qalpha1(1)-mean(alphasim1), qalpha1(2)-mean(alphasim1), ...
         'horizontal','LineWidth',1,'color','b', 'HandleVisibility', 'off')
     


