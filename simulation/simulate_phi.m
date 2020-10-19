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

    %ubroadening is the spectral resolution 
    %of the measurements divided by bin width u of the spectra.
    ubroadening = 3; 

    %L regularization matrix
    L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

    %Observation angles
    phi = [60 80];
    thirdangle = [0:5:90];
    
    %Set the third angle.
    phi(3) = thirdangle(idx);

    [b, binfo] = biMaxb(ustruct,phi);
    [b_noisy, ~, e] = add_noise(b,0.01);

    A = biMaxA(ubroadening, xinfo, binfo);

    welford.keep = 1;           %No trimming
    
    %welford.nburnin = 1;
    %nsim = 10;
    welford.nburnin = 500;     %.5k burnin
    nsim = 5000;                %5k samples

    [x, alpha, delta, lambda, info] = NNHGS(A,b_noisy,L,nsim,welford);

    %Save all of the current variables to a .mat file
    foldername = angle_to_string(phi);
    save(strcat('./sim_angles2/', angle_to_string(phi), '.mat'))

    %Collect the angles in a string
    function s = angle_to_string(phi)
        s = sprintf('phi_%d_%d_%d',phi(1),phi(2),phi(3));
    end
end
