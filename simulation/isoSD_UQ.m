clear, clc, close all

%Set the simulation setup.
simname = 'isoSD_UQ';
nsim    = 5000; %Number of samples simulated per run.
disp_waitbar = 1;

if ~isfolder(simname)
    mkdir(simname)
else
    error('Folder already exists')
end

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate on this grid
[F, isoSDFinfo] = isoSDF(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, isoSDSinfo] = isoSDS(ustruct,phi);
%S_noisy = add_noise(S,0.01); 
%Some of the elements in S are negative, and therefore the noise will be
%complex.
S_noisy = S + 0.05*randn(size(S));

S = S_noisy;
W = isoSDW(3,isoSDFinfo,isoSDSinfo);

save(strcat('./',simname,'/setup.mat'))

%% UQ on isoSD case.


while 1
   clearvars -except W S alpha nsim simname disp_waitbar
   alpha = 2.2230e+09; %alpha0 is pretty important and dependent on the noise.
   [F_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(W,S,alpha,nsim,disp_waitbar);
   
   save(strcat('./',simname,'/sim',num2str(nextsimnum(simname)),'.mat'), ...
        'F_sim','del_sim','lam_sim','alph_sim')
end



