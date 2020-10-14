%% biMax Angle Simulation 14-10-2020
clear, clc, close all

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();
vparadim = length(gridinfo.vpara_ax);
vperpdim = length(gridinfo.vperp_ax);

%Evaluate the bi-Maxwellian on this grid with default values.
[x_true, xinfo] = biMaxx(vpara, vperp); x_true = x_true(:);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;

%ubroadening is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
ubroadening = 3; 

%L regularization matrix
L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

%Observation angles
phi = [60 80];
thirdangle = [0:5:90];

for i=1:length(thirdangle)
    tic
    phi(3) = thirdangle(i);
    disp(angle_to_string(phi))
    
    [b, binfo] = biMaxb(ustruct,phi);
    [b_noisy, ~, e] = add_noise(b,0.01);
    
    A = biMaxA(ubroadening, xinfo, binfo);
    
    welford.keep = 1;           %No trimming
    welford.nburnin = 1000;     %1k burnin
    nsim = 10000;               %10k samples
    
    [x, alpha, delta, lambda, info] = NNHGS(A,b_noisy,L,nsim,welford);
    
    %Save all of the current variables to a .mat file
    foldername = angle_to_string(phi);
    mkdir(strcat('./sim_angles/',foldername))
    save(strcat('./sim_angles/',foldername,'/sim.mat'))
    toc
end

%Collect the angles in a string
function s = angle_to_string(phi)
    s = sprintf('angle_%d_%d_%d',phi(1),phi(2),phi(3));
end