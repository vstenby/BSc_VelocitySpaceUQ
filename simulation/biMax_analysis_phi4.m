clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_phi4')
    system('scp -r s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/biMax_sim_phi4 biMax_sim_phi4');
end

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('biMax_sim_phi4')

files = {dir('biMax_sim_phi4/*.mat').name}; 

phi_analysis_table(files)
phi_analysis_plot(files)