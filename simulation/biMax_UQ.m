%% Drifting bi-Maxwellian UQ
%Commented out so it doesn't reset on HPC.
clear, clc, close all

%Set the simulation setup.
simname = 'biMax_UQ';
nsim    = 5000; %Number of samples simulated per run.
disp_waitbar = 1;

if ~isfolder(simname)
    mkdir(simname)
elseif input('Should simulations be added to the same folder? ') == 1
    %Nothing should happen.
else
    error('Folder already exists')
end

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate the bi-Maxwellian on this grid with default values.
[F, biMaxFinfo] = biMaxF(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, biMaxSinfo] = biMaxS(ustruct,phi);

[S_noisy, s] = add_noise(S,0.01);

S = S_noisy;

W = biMaxW(3,biMaxFinfo,biMaxSinfo);

save(strcat('./',simname,'/setup.mat'))

%% Do the sampling
while 1
   clearvars -except W S alpha nsim simname disp_waitbar
   alpha = 2.8e8; %alpha0
   
   [F_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(W,S,alpha,nsim,disp_waitbar);
   
   save(strcat('./',simname,'/sim',num2str(nextsimnum(simname)),'.mat'), ...
        'F_sim','del_sim','lam_sim','alph_sim')
end
