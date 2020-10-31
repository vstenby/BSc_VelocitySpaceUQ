clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_rhs')
    HPCDownload('biMax_sim_rhs','biMax_sim_rhs','s174483');
end

addpath('biMax_sim_rhs')
files = {dir('biMax_sim_rhs/*.mat').name}; 

rhs_analysis_plot(files)