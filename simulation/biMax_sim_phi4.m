%% Bi-Maxwellian Distribution Simulation (Simulate the 4th angle)
function biMax_sim_phi4(sim_idx)

    %Add dependencies to path.
    addpath(genpath('../functions'))
    addpath(genpath('../../aux'))
    
    % - - - - - - - - - - Simulation setup  - - - - - - - - - - - - - - - -
    % In this simulation, the first two observation angles are 60 and 80,
    % and the last one is in the range 1 to 90.
    
    %Observation angles.
    phi  = [60 80 20];
    phi4 = [5:5:90]; phi4 = [1 phi4];
    phi(4) = phi4(sim_idx); clear phi4
    
    nsim = 5000;
    nburnin = 500;
  
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    rhs = 2;
    
    %Parameters for the (vpara,vperp)-grid.
    vparamin=-4e6;
    vparamax=4e6;
    vperpmin=1e4;
    vperpmax=4e6;

    %Construct the u-vector.
    u = construct_uvec('umax', 5*1e6, 'du', 1e5);

    %Construct the normal grid and our true solution.
    vparadim = 35; vperpdim = 17;
    [vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
    xtrue = biMaxx(vpara,vperp);

    %Forward model matrix A
    A = transferMatrix(vpara,vperp,phi,u);

    switch rhs
        case 1
            [b, e] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 
        case 2
            vparadimfine = 100; vperpdimfine = 50;
            [vparafine, vperpfine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
            xfine = biMaxx(vparafine, vperpfine); 
            Afine = transferMatrix(vparafine,vperpfine,phi,u); 
            [b, e] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.
        case 3
            gu = biMaxg(u,phi);
            du = u(2)-u(1);
            b = du*gu;
            [b, e] = add_noise(b);
    end

    %Normalize with this noise.
    [A, b] = error_normalization(A,b,e);    

    %Regularization matrix L (1st order Tikhonov)
    L = reguL(vperpdim, vparadim);

    %A vector of alphas for finding the optimal regularization parameter.
    alpha_relerr = logspace(-9,-4,100);

    %Find the relative error for 0th order and 1st order.
    [~, ~, r0th] = TikhNN(A, b, alpha_relerr, [], 'return_relerr', true, 'x_true', xtrue);
    [~, ~, r1st] = TikhNN(A, b, alpha_relerr,  L, 'return_relerr', true, 'x_true', xtrue);

    %Find the alpha with the smallest relative error.
    [minr0,idx0] = min(r0th); optalpha_0th = alpha_relerr(idx0);
    [minr1,idx1] = min(r1st); optalpha_1st = alpha_relerr(idx1);

    %Find the corresponding solution.
    xopt0 = TikhNN(A,b,optalpha_0th);
    xopt1 = TikhNN(A,b,optalpha_1st,L);

    %Do the sampling.
    [xNNHGS0, alphasim0, deltasim0, lambdasim0, NNHGS0info] = NNHGS(A,b,[],nsim,'solver', 'lsqnonneg', 'welford',true,'nburnin',nburnin);
    [xNNHGS1, alphasim1, deltasim1, lambdasim1, NNHGS1info] = NNHGS(A,b,L,nsim,'solver','lsqnonneg','welford',true,'nburnin',nburnin);

    %Calculate reconstruction with mean(alpha)
    xsamplealpha0 = TikhNN(A,b,mean(alphasim0));
    xsamplealpha1 = TikhNN(A,b,mean(alphasim1),L);

    %Empiric confidence intervals for alpha.
    CIalpha0 = quantile(alphasim0,[0.025, 0.975]);
    CIalpha1 = quantile(alphasim1,[0.025, 0.975]);

    %Calculate the caxis for all solutions and standard deviation of solution.
    caxis_mu  = [min([xtrue(:) ; xNNHGS0(:,1) ; xNNHGS1(:,1)]), ...
                 max([xtrue(:) ; xNNHGS0(:,1) ; xNNHGS1(:,1)])];

    caxis_std = [min([xNNHGS0(:,2) ; xNNHGS1(:,2)]), ...
                 max([xNNHGS0(:,2) ; xNNHGS1(:,2)])];
             
    outputpath = sprintf('./biMax_sim_phi4/phi4_60_80_20_%02d.mat',phi(4));
    save(outputpath)
end