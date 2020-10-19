%% Script made for analyzing the angles.
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%% Analyze the different angles.
addpath('sim_angles')

%Observation angles
phi = [60 80];
thirdangle_load = [0:5:90];

for i=1:length(thirdangle_load)
    phi(3) = thirdangle_load(i);
    loadpath = strcat('./sim_angles/',angle_to_string(phi),'/sim.mat');
    load(loadpath);
    geweketest([alpha,delta,lambda])
    xload(:,1:2,i) = x;
end
disp('Done loading')

for i=1:size(xload,3)
    phi(3) = thirdangle_load(i);
    figure(i)
    sgtitle(angle_to_string(phi))
    subplot(1,2,1)
    showDistribution(xload(:,1,i),gridinfo);
    %imagesc(x(:,1,i)); axis xy; caxis(caxis_mu);
    subplot(1,2,2)
    showDistribution(xload(:,2,i),gridinfo);
    %imagesc(x(:,2,i)); axis xy; caxis(caxis_std);
end

%%


%Collect the angles in a string
function s = angle_to_string(phi)
    s = sprintf('angle_%d_%d_%d',phi(1),phi(2),phi(3));
end



