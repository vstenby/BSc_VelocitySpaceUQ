%% biMax Angle Simulation 14-10-2020
function simulate_phi(idx)
    %Add dependencies.
    addpath(genpath('../functions'))
    addpath(genpath('../../aux'))

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

    %L regularization matrix
    %0th order
    %L = []; 
    %1st order
    L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.
    
    %Observation angles
    phi = [60 80 20];
    fourthangle = [0:5:90];
    
    %Set the third angle.
    phi(4) = fourthangle(idx);
    
    [b, binfo] = biMaxb(ustruct,phi);
    [b_noisy, ~, e] = add_noise(b,0.01);
    
    u = construct_uvec(ustruct);
    A = transferMatrix(vpara,vperp,phi,u);

    welford.nburnin = 500;     %.5k burnin
    nsim = 5000;                %5k samples

    
    %Save all of the current variables to a .mat file
    foldername = angle_to_string(phi);
    save(strcat('./sim_4angles_1st/', angle_to_string(phi), '.mat'))
    
    %Collect the angles in a string
    function s = angle_to_string(phi)
        s = sprintf('phi_%d_%d_%d_%d',phi(1),phi(2),phi(3),phi(4));
    end
end
