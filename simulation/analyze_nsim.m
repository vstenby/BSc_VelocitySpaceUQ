%% Script for analyzing the nsim numerical experiment.
clear, clc, close all

%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('./nsim')

load('./nsim/nangle2.mat');

%% Show the chains
clc

%For 2 angles, 90 samples is enough. 
k = 100;
geweketest([info.alpha(1:k), info.delta(1:k), info.lambda(1:k)], 10)



