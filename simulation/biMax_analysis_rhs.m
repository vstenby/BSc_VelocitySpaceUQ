clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_rhs')
    system('scp -r s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/biMax_sim_rhs biMax_sim_rhs');
end

addpath('biMax_sim_rhs')
files = {dir('biMax_sim_rhs/*.mat').name}; 

rhs_analysis_plot(files)