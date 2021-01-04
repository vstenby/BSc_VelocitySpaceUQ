%% Tikhonov Reconstruction with the Bi-Maxwellian Velocity Distribution
clear, clc, close all
% ------------------- Settings for the example -------------------------- %

% Set the right hand side as either 'crime', 'no crime' or 'analytic'.
rhs = "analytic"; 

% nalpha is the number of alphas that TikhNN tests for both the 0th and
% 1st order Tikhonov. Therefore, nalpha = 50 would mean we evaluate 100
% alphas. Consider setting to a lower number if it takes too long.
nalpha = 100;

% ----------------------------------------------------------------------- %

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

[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim, ...
'vperpmin',vperpmin,'vperpmax',vperpmax, ... 
'vparamin',vparamin,'vparamax',vparamax);

%Generate the true solution.
x_true = biMaxx(vpara,vperp);

%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

%Generate the right hand side
switch rhs
    case 'crime'
        %With crime, we use A both for the forward and backwards on the
        %same grid. generate_noisy_b generates a noisy measurement as
        %b = A*x + e
        [b, e] = generate_noisy_b(A,x_true);
    case 'no crime'
        %With no crime, we use a finer grid for the forward, and then
        %coarser grid backwards. 
        vparadimfine = 100; vperpdimfine = 50;
        [vparafine, vperpfine, ginfofine] = construct_vgrid(vparadimfine, ...
        vperpdimfine,'vparamax',vparamax,'vparamin', ...
        vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
        xfine = biMaxx(vparafine, vperpfine); 
        Afine = transferMatrix(vparafine,vperpfine,phi,u); 
        [b, e] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.
    case 'analytic'
        %Analytic right hand side
        gu = biMaxg(u,phi);
        du = u(2)-u(1);
        b = du*gu;
        [b, e] = add_noise(b);
    otherwise
        error('Wrong right hand side specified.')
end

%Here, we normalize with the error. 
[A, b] = error_normalization(A,b,e);

%Regularization matrix L (1st order Tikhonov)
L = reguL(vpara, vperp);

%0th order Tikhonov values
alphas0 = logspace(-15,-8,nalpha);
alphas1 = logspace(-4,4,nalpha);

%Find the relative error for 0th order and 1st order.
[x0, ~, r0] = TikhNN(A, b, alphas0, [], 'return_relerr', true, 'x_true', x_true);
[x1, ~, r1] = TikhNN(A, b, alphas1,  L, 'return_relerr', true, 'x_true', x_true);

%Fetch information about the best solution (relative error-wise)
[r0min, idx0] = min(r0); alpha0opt = alphas0(idx0); x0opt = x0(:,idx0);
[r1min, idx1] = min(r1); alpha1opt = alphas1(idx1); x1opt = x1(:,idx1);

%Plots
figure
sgtitle(sprintf('Optimal solution for 0th order Tikhonov, r = %.2f', r0min))
subplot(2,2,[1,3])
semilogx(alphas0, r0, 'k-'); xlabel('\alpha'); ylabel('Relative error')
hold on; plot(alpha0opt, r0min, 'r.', 'MarkerSize', 15)
subplot(2,2,2)
showDistribution(x_true, ginfo); title('True solution')
subplot(2,2,4)
showDistribution(x0opt, ginfo); title('Tikhonov solution')

figure
sgtitle(sprintf('Optimal solution for 1st order Tikhonov, r = %.2f', r1min))
subplot(2,2,[1,3])
semilogx(alphas1, r1, 'k-'); xlabel('\alpha'); ylabel('Relative error')
hold on; plot(alpha1opt, r1min, 'r.', 'MarkerSize', 15)
subplot(2,2,2)
showDistribution(x_true, ginfo); title('True solution')
subplot(2,2,4)
showDistribution(x1opt, ginfo); title('Tikhonov solution')

figure
subplot(1,3,1)
showDistribution(x_true,ginfo); title('True solution')
subplot(1,3,2)
showDistribution(x0opt,ginfo); title('0th order Tikhonov')
subplot(1,3,3)
showDistribution(x1opt,ginfo); title('1st order Tikhonov')