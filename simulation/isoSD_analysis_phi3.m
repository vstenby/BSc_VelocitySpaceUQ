clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('isoSD_sim_phi3')
    system('scp -r s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/isoSD_sim_phi3 isoSD_sim_phi3');
end

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('isoSD_sim_phi3')

files = {dir('isoSD_sim_phi3/*.mat').name}; 

phi_analysis_table(files)
phi_analysis_plot(files)