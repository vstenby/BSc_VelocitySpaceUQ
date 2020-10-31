clear, clc, close all

HPCDownload('biMax_sim_nsim','biMax_sim_nsim','s174483')

%%

files = {dir('biMax_sim_nsim/*.mat').name}; 
nfiles = length(files);
addpath('biMax_sim_nsim')

nsim_analysis_table(files)
