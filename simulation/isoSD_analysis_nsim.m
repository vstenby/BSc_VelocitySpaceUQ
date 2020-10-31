clear, clc, close all

HPCDownload('isoSD_sim_nsim','isoSD_sim_nsim','s174483')

%%

files = {dir('isoSD_sim_nsim/*.mat').name}; 
nfiles = length(files);
addpath('isoSD_sim_nsim')

nsim_analysis_table(files)
