function data_sim(sim_idx)
% Simulates the different setups.
%clear, clc, close all

%sim_idx = 1;
% We have 7 different timestamps, and for each timestamp we have 4 rhs.
% This means a total of 28 simulations.

% -- Make all combinations of simulations.
sims = combvec(1:7, 1:4);

sim = sims(:,sim_idx); ts = sim(1); rhs = sim(2);

loadfile = sprintf('./data/ts%d.mat',ts);
outfile = sprintf('./data_sim/ts_%02d_rhs_%02d.mat',ts,rhs);

switch rhs
    case 1
    %Inverse crime
    load(loadfile, 'f_TR_coarse', 'L1p', 'L1E', 'L0b', 'transf_matrix', 'S_TR_invcrime_n');
    A = transf_matrix;
    b = double(S_TR_invcrime_n);
    
    case 2
    %Non-inverse crime    
    load(loadfile, 'f_TR_coarse', 'L1p', 'L1E', 'L0b', 'transf_matrix', 'S_TR_n');
    A = transf_matrix;
    b = double(S_TR_n);
    
    case 3
    %FIDA sim
    load(loadfile, 'f_TR_coarse', 'L1p', 'L1E', 'L0b', 'transf_matrix', 'fida_total_n');
    A = transf_matrix;
    b = double(fida_total_n);
    
    case 4
    %Real data
    load(loadfile, 'f_TR_coarse', 'L1p', 'L1E', 'L0b', 'transf_matrix_broad', 'S_data', 'err_data');
    A = transf_matrix_broad;
    b = double(S_data);
    
    %Normalize with the error.
    [A, b] = error_normalization(A, b, err_data);
end

%True solution
xtrue = f_TR_coarse;

%Save L into a struct.
L.L1p = L1p; 
L.L1E = L1E;
L.L0b = L0b;

% -- Tikhonov for different values of alpha --
alpha = logspace(-6,5,20); 
[xalpha, ~, relerr] = TikhNN(A, b, alpha, L, 'return_relerr', true, 'x_true', xtrue);
% -- End of Tikhonov --

nsim = 5000;
nburnin = 500;

% -- Nonnegative Gibbs Sampler --
[xsim, alphasim, deltasim, lambdasim, info] = NNHGS(A, b, L, nsim, 'welfordsalgorithm', true, 'nburnin', nburnin);

save(outfile)
end