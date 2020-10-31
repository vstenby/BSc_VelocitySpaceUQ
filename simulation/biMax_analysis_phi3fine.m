clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_phi3fine')
    system('scp s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/biMax_sim_phi3fine/*.mat biMax_sim_phi3fine');
end

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('biMax_sim_phi3fine')

files = {dir('biMax_sim_phi3fine/*.mat').name}; 

phi_analysis_table(files)
phi_analysis_plot(files)