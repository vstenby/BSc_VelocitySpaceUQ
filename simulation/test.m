clear, clc, close all
%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

folder = 'sim_angles_1st';
%sim1 = load_simulation(folder);

sim1 = UQSim(folder, 'drop', {'phi_60_80_45.mat', 'phi_60_80_80.mat', 'phi_60_80_60.mat', 'phi_60_80_90.mat'}, 'keep1st', {'gridinfo', 'x_true'});

%[x, x_true, alpha, delta, lambda, gridinfo, filename]  = load_simulation('sim_angles_0th');

sim1.slideshow()