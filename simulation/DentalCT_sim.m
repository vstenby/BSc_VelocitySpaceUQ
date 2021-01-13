function DentalCT_sim(sim_idx)

    addpath(genpath('../other/DentalCT'))
    addpath(genpath('../functions'))
    addpath(genpath('../../aux'))

    thetas = [180-95:2.5:180-95/2; ...
              180-95:5:180; ...
              180-95:10:180+95; ...
              180-95/2:5:180+95/2; ...
              180:5:180+95; ...
              180+95/2:2.5:180+95];

    theta = thetas(sim_idx,:); clear thetas
    
    nsim = 5500;
    
    % --- Set up the problem ---
    N = 30; 
    p = 35;
    R = 6;
    dw = 1.5;
    sd = 9;

    A = fanlineartomo(N,theta,p,R,dw,sd);
    x = dental_phantom(N); 
    Ax = A*x(:);
    SNR = 50;
    err_lvl = 100/SNR;
    noise = err_lvl/100 * norm(Ax(:)) / sqrt(size(A,1));

    b = Ax + noise*randn(size(Ax));
    % --- Go on to Tikhonov ---

    I = @(d) speye(d);
    D = @(d) spdiags([zeros(d-1,1) -ones(d-1,1) ones(d-1,1)],-1:1,d-1,d);

    L = [kron(D(30), I(30)) ; kron(I(30),D(30))];

    alphavec = logspace(-4,1,100);
    [~,~,r0] = TikhNN(A,b,alphavec,[],'return_relerr',true,'x_true',x(:));
    [~,~,r1] = TikhNN(A,b,alphavec,L,'return_relerr',true,'x_true',x(:));

    % -- Go onto Gibbs Sampler --

    %Find the alpha with the smallest relative error.
    [minr0,idx0] = min(r0); optalpha_0th = alphavec(idx0);
    [minr1,idx1] = min(r1); optalpha_1st = alphavec(idx1);

    %Find the corresponding solution.
    xopt0 = TikhNN(A,b,optalpha_0th);
    xopt1 = TikhNN(A,b,optalpha_1st,L);

    %Do the sampling.
    [xNNHGS0, alphasim0, deltasim0, lambdasim0, NNHGS0info] = NNHGS(A,b,[],nsim);
    [xNNHGS1, alphasim1, deltasim1, lambdasim1, NNHGS1info] = NNHGS(A,b,L,nsim);

    %Calculate reconstruction with mean(alpha)
    xsamplealpha0 = TikhNN(A,b,mean(alphasim0));
    xsamplealpha1 = TikhNN(A,b,mean(alphasim1),L);

    %Empiric confidence intervals for alpha.
    qalpha0 = quantile(alphasim0,[0.025, 0.975]);
    qalpha1 = quantile(alphasim1,[0.025, 0.975]);

    outputpath = sprintf('./DentalCT_sim/DentalCT_setup%02d.mat',sim_idx);
    save(outputpath) 

end