clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('isoSD_sim_rhs')
    system('scp -r s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/isoSD_sim_rhs isoSD_sim_rhs');
end

addpath('isoSD_sim_rhs')
files = {dir('isoSD_sim_rhs/*.mat').name}; 

rhs_analysis_plot(files)