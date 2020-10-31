clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_phi3')
    HPCDownload('biMax_sim_phi3','biMax_sim_phi3','s174483');
end

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('biMax_sim_phi3')

files = {dir('biMax_sim_phi3/*.mat').name}; 

phi_analysis_table(files)
phi_analysis_plot(files)

