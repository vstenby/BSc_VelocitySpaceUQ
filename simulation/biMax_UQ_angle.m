%%  Bi-Maxwellian UQ of different angles.
%
%   This is the bi-Maxwellian where we do UQ for different angles.
%

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Set the simulation setup.
simname = 'biMaxUQ_angles'; %Foldername.
nsim    = 2000;

phis = {[40, 70], [40, 85], ... 
        [70, 85], ...
        [10, 20, 40], [10, 20, 70], [10, 20, 85], ...
        [10, 40, 70], [10, 40, 85], ...
        [20, 40, 70], [20, 40, 85], ...
        [20, 70, 85], ...
        [40, 70, 85], ...
        [10, 20, 40, 70], [10, 20, 70, 85]};
    
if ~isfolder(simname)
    mkdir(simname)
elseif input('Should simulations be added to the same folder? ') == 1
    %Nothing should happen.
else
    error('Folder already exists')
end

for i=1:length(phis)
    clearvars -except nsim phis simname i
    %Observation angles.
    phi = phis{i};
    
    %Subfoldername
    subfoldername = strcat('./', simname, '/', num2str(phi));
    mkdir(subfoldername)
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

    [b, binfo] = biMaxb(ustruct,phi);

    b_noisy = add_noise(b,0.01);

    ubroadening = 3;
    A = biMaxA(ubroadening,xinfo, binfo);

    %This is for the 1st order Tikhonov.
    L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

    %L has to be square for UQ.
    %L = chol(L'*L);
    
    save(strcat(subfoldername, '/setup.mat'))
    
    % 3a) UQ for the 0th order Tikhonov formulation
    [xUQ0th, del0th, lam0th, alph0th] = NNHGS_UQ(A, b_noisy, 0, [], nsim, 1);
    save(strcat(subfoldername,'/0thTikhUQ.mat'), 'xUQ0th','del0th','lam0th','alph0th')
    
    % 3b) UQ for the 1st order Tikhonov formulation
    [xUQ1st, del1st, lam1st, alph1st] = NNHGS_UQ(A, b_noisy, 0, L, nsim, 1);
    save(strcat(subfoldername,'/1stTikhUQ.mat'), 'xUQ1st','del1st','lam1st','alph1st')
end
    