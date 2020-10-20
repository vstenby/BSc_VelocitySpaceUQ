%% Analysis of the "simulate_phi" experiment.

%Load the data.
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

foldername = 'sim_angles_1st'; addpath(foldername);

files = dir(strcat(foldername,'/*.mat'));
files = {files.name};

files = sortfiles(files);

load(strcat(foldername,'/',files{1}),'x_true', 'phi', 'alpha', 'phi');
nsimulations = length(files);
nsamps = size(alpha,1);


phi_sim    = zeros(nsimulations,length(phi));
pgeweke    = zeros(nsimulations,1);
alpha_sim  = zeros(nsamps,nsimulations);
delta_sim  = zeros(nsamps,nsimulations);
lambda_sim = zeros(nsamps,nsimulations);
x_sim      = zeros(size(x_true,1), 2, nsimulations);
phi_sim    = zeros(nsimulations, length(phi));

for i=1:length(files)
    load(strcat(foldername,'/',files{i}), 'alpha', 'delta', 'lambda','x','phi','gridinfo')
    pgeweke(i) = geweketest([alpha delta lambda]);
    alpha_sim(:,i) = alpha;
    delta_sim(:,i) = delta;
    lambda_sim(:,i) = lambda;
    x_sim(:,:,i) = x;
    phi_sim(i,:) = phi;
end
disp('Done loading')

%% Show the different reconstructions.
clim_mu  = [min(reshape(x_sim(:,1,:),1,[])) ...
           max(reshape(x_sim(:,1,:),1,[]))];
       
clim_std = [min(reshape(x_sim(:,2,:),1,[])) ...
           max(reshape(x_sim(:,2,:),1,[]))];

close all
figure('units','normalized','outerposition',[0 0 1 1])
idx = 1;

while idx ~= 0 && idx ~= nsimulations+1
    subplot(1,2,1)
    showDistribution(x_sim(:,1,idx),gridinfo, clim_mu)
    title(strcat(phi2str(phi_sim(idx,:)),' mean'))
    subplot(1,2,2)
    showDistribution(x_sim(:,2,idx),gridinfo, clim_std)
    title(strcat(phi2str(phi_sim(idx,:)),' std'))
    drawnow()
    k = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    switch value
        case 28 %Left arrow
            idx = idx - 1;
        case 29 %Right arrow
            idx = idx + 1;
        otherwise
            close all
            break
    end
end

%% Plot all of the svd values
close all
idx = 1;
figure('units','normalized','outerposition',[0 0 1 1])
while idx ~= 0 && idx ~= nsimulations+1
    load(strcat('./',foldername,'/',files{idx}),'A','phi');
    [~,S] = svd(A);
    semilogy(diag(S)); title(phi2str(phi));
    drawnow()
    hold off
    k = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    switch value
        case 28 %Left arrow
            idx = idx - 1;
        case 29 %Right arrow
            idx = idx + 1;
        otherwise
            close all
            break
    end
end

%% Plot the alpha parameters vs. relative error.
alphavec = logspace(7, 11, 20);
close all
figure('units','normalized','outerposition',[0 0 1 1])
idx = 1;
while idx ~= 0 && idx ~= nsimulations+1
    load(strcat('./',foldername,'/',files{idx}))
    x1st = GPCG_TikhNN(A,b_noisy,alphavec,L);
    [r1, idx1] = relerr(x_true, x1st);
    qalpha = quantile(alpha, [0, 0.25, 0.5, 0.75, 1]);
    %Make the plot
    semilogx(alphavec,r1, 'g--')
    hold on
    plot(alphavec(idx1),r1(idx1), 'g.', 'MarkerSize',15)
    hold on
    
    %Draw the quantiles
    xline(qalpha(1), 'r--');
    xline(qalpha(2), 'r-');
    xline(qalpha(3), 'k-');
    xline(qalpha(4), 'b-');
    xline(qalpha(5), 'b--');
   
    drawnow()
    hold off
    k = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    switch value
        case 28 %Left arrow
            idx = idx - 1;
        case 29 %Right arrow
            idx = idx + 1;
        otherwise
            close all
            break
    end
end

%% Auxil. functions

function n = extract_n(filename)
idx1 = regexp(filename,'_'); idx1 = idx1(end);
idx2 = regexp(filename,'.mat'); 
n = str2num(filename(idx1+1:idx2-1));
end

function files_sorted = sortfiles(files)
% sort the files by phi.
ns = zeros(1,length(files));
for i=1:length(files)
   ns(i) = extract_n(files{i});
end
[~, idx] = sort(ns);
files_sorted = {files{idx}};
end

function s = phi2str(phi)
s = sprintf('%02d %02d %02d', phi(1), phi(2), phi(3));
end


