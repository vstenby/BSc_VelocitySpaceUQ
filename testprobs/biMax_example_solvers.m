%% Bi-Maxwellian Distribution example.
clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[10 20 40 70 85];

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

gu = biMaxg(u,phi);
du = u(2)-u(1);
b = du*gu;
[b, e] = add_noise(b);

%Normalize with this noise.
[A, b] = error_normalization(A,b,e);    

%Regularization matrix L (1st order Tikhonov)
L = reguL(vperpdim, vparadim);

%A vector of alphas for finding the optimal regularization parameter.
alpha_relerr = logspace(-9,-4,20);

if ~check_mosek()
    addpath(genpath('../../aux'))
    if ~check_mosek(), error('mosek not loaded'), end
end
    
fprintf('Starting 0th order mosek lsqlin, ')
tic
[~, ~, r0_lsqlin_mosek] = TikhNN(A, b, alpha_relerr, [], 'return_relerr', true, 'x_true', xtrue);
fprintf('done in %f seconds.\n',toc)

fprintf('Starting 1st order mosek lsqlin, ')
tic
[~, ~, r1_lsqlin_mosek] = TikhNN(A, b, alpha_relerr,  L, 'return_relerr', true, 'x_true', xtrue);
fprintf('done in %f seconds.\n',toc)

fprintf('Starting 0th order mosek lsqnonneg, ')
tic
[~, ~, r0_lsqnonneg_mosek] = TikhNN(A, b, alpha_relerr, [], 'return_relerr', true, 'x_true', xtrue, 'solver', 'lsqnonneg');
fprintf('done in %f seconds.\n',toc)

fprintf('Starting 1st order mosek lsqnonneg, ')
tic
[~, ~, r1_lsqnonneg_mosek] = TikhNN(A, b, alpha_relerr,  L, 'return_relerr', true, 'x_true', xtrue, 'solver', 'lsqnonneg');
fprintf('done in %f seconds.\n',toc)


rmpath(genpath('../../aux')); if check_mosek(), error('Unavailable to remove mosek'), end

fprintf('Starting 0th order MATLAB lsqlin, ')
tic
[~, ~, r0_lsqlin_matlab] = TikhNN(A, b, alpha_relerr, [], 'return_relerr', true, 'x_true', xtrue);
fprintf('done in %f seconds.\n',toc)

fprintf('Starting 1st order MATLAB lsqlin, ')
tic
[~, ~, r1_lsqlin_matlab] = TikhNN(A, b, alpha_relerr,  L, 'return_relerr', true, 'x_true', xtrue);
fprintf('done in %f seconds.\n',toc)

%% Plot the relative errors
semilogx(r0_lsqlin_matlab, 'b-')
hold on
semilogx(r1_lsqlin_matlab, 'b-.')
semilogx(r0_lsqlin_mosek,'r--')
semilogx(r1_lsqlin_mosek,'r--')
semilogx(r0_lsqnonneg_mosek,'g.')
semilogx(r1_lsqnonneg_mosek,'g.')
