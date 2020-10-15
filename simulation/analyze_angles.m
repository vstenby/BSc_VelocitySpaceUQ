%% Script made for analyzing the angles.
clear, clc, close all

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%% Analyze the different angles.
addpath('sim_angles')
load('./sim_angles/angle_60_80_5/sim.mat');

geweketest([alpha,delta,lambda])

figure
subplot(1,2,1)
showDistribution(x(:,1),gridinfo)
subplot(1,2,2)
showDistribution(x(:,2),gridinfo)



