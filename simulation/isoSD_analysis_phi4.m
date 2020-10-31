clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('isoSD_sim_phi4')
    system('scp -r s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/isoSD_sim_phi4 isoSD_sim_phi4');
end

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('isoSD_sim_phi4')

%files = {dir('isoSD_sim_phi4/*.mat').name}; 
files = {'phi4_60_80_20_90.mat'};
phi_analysis_table(files)
phi_analysis_plot(files)

