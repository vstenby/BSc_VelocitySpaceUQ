%% Isotropic Slowing Down Example
%
% In this example, three different right hand sides are constructed.
% b1 : Inverse crime (i.e. coarse grid forward)
% b2 : Non-inverse crime (i.e. fine grid forward)
% b3 : Analytic rhs.
%
% After we have seen that the three right hands are indeed very similar,
% then we do the Tikhonov Regularisation for various alpha values for the
% three right hand sides. 

clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Parameters for the (vpara,vperp)-grid.
vparamin = -1.4e7;
vparamax = 1.4e7;
vperpmin = 1e5;
vperpmax = 1.4e7;

%Observation angles.
phi = [10 20 40 70 85];

%Construct the u-vector.
%u = construct_uvec('Emin', 10e3, 'Emax', 4e6, 'udim', 100);
u = construct_uvec('Emax', 4e6, 'udim', 100);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = isoSDx(vpara,vperp);

%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

%Generate the inverse crime 
[b1, e1] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 

%Non-inverse crime
vparadimfine = 100; vperpdimfine = 50;
[vparafine, vperpfine, ginfofine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
xfine = isoSDx(vparafine, vperpfine); 
Afine = transferMatrix(vparafine,vperpfine,phi,u); 
[b2, e2] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.

%Analytic right hand side
gu = isoSDg(u,phi);
du = u(2)-u(1);
b3 = du*gu;
[b3, e3] = add_noise(b3);

%Regularization matrix L (1st order Tikhonov)
L = reguL(vpara, vperp);

%Now we wish to do the Tikhonov reconstruction for each of the three cases.
%Here, we correct A and b with the error. 

[A1, b1] = error_normalization(A,b1,e1);  
[A2, b2] = error_normalization(A,b2,e2);
[A3, b3] = error_normalization(A,b3,e3);

%% -- 0th order Tikhonov --
disp('Solving 0th order Tikhonov problems...')
%A vector of alphas for finding the optimal regularization parameter.
alpha_relerr_0th = logspace(-12,-5,100);

%Find the relative error for 0th order and 1st order.
[~, ~, r0th_crime]    = TikhNN(A1, b1, alpha_relerr_0th, [], 'return_relerr', true, 'x_true', xtrue);
[~, ~, r0th_nocrime]  = TikhNN(A2, b2, alpha_relerr_0th, [], 'return_relerr', true, 'x_true', xtrue);
[~, ~, r0th_analytic] = TikhNN(A3, b3, alpha_relerr_0th, [], 'return_relerr', true, 'x_true', xtrue);

semilogx(alpha_relerr_0th, r0th_crime);
hold on
semilogx(alpha_relerr_0th, r0th_nocrime);
hold on
semilogx(alpha_relerr_0th, r0th_analytic);
legend('Inverse crime', 'No crime', 'Analytic', 'location', 'southwest', 'FontSize',15)
xlabel('$\alpha$','Interpreter','Latex','FontSize',15)
ylabel('Relative error','FontSize',15)

%Find the corresponding alphas with the smallest relative errors.
[minr0_crime,idx0_crime]       = min(r0th_crime);     optalpha_0th_crime    = alpha_relerr_0th(idx0_crime);
[minr0_nocrime,idx0_nocrime]   = min(r0th_nocrime);   optalpha_0th_nocrime  = alpha_relerr_0th(idx0_nocrime);
[minr0_analytic,idx0_analytic] = min(r0th_analytic);  optalpha_0th_analytic = alpha_relerr_0th(idx0_analytic);

%Find the corresponding solutions
xopt0_crime    = TikhNN(A1,b1,optalpha_0th_crime);
xopt0_nocrime  = TikhNN(A2,b2,optalpha_0th_nocrime);
xopt0_analytic = TikhNN(A3,b3,optalpha_0th_analytic);

%%
figure
semilogx(alpha_relerr_0th, r0th_crime);
hold on
semilogx(alpha_relerr_0th, r0th_nocrime);
hold on
semilogx(alpha_relerr_0th, r0th_analytic);
legend('Inverse crime', 'No crime', 'Analytic', 'location', 'southwest', 'FontSize',15)
xlabel('$\alpha$','Interpreter','Latex','FontSize',15)
ylabel('Relative error','FontSize',15)

figure
showDistribution(xopt0_crime,ginfo); title(sprintf('Inverse crime, relative error: %.4f', minr0_crime))

figure
showDistribution(xopt0_nocrime,ginfo); title(sprintf('No crime, relative error: %.4f', minr0_nocrime))

figure
showDistribution(xopt0_analytic,ginfo); title(sprintf('Analytic, relative error: %.4f', minr0_analytic))

%% -- 1st order Tikhonov -- 
disp('Solving 1st order Tikhonov problems...')
%A vector of alphas for finding the optimal regularization parameter.
alpha_relerr_1st = logspace(1,7,100);

%Find the relative error for 0th order and 1st order.
[~, ~, r1st_crime]    = TikhNN(A1, b1, alpha_relerr_1st, L, 'return_relerr', true, 'x_true', xtrue);
[~, ~, r1st_nocrime]  = TikhNN(A2, b2, alpha_relerr_1st, L, 'return_relerr', true, 'x_true', xtrue);
[~, ~, r1st_analytic] = TikhNN(A3, b3, alpha_relerr_1st, L, 'return_relerr', true, 'x_true', xtrue);

figure
semilogx(alpha_relerr_1st, r1st_crime);
hold on
semilogx(alpha_relerr_1st, r1st_nocrime);
hold on
semilogx(alpha_relerr_1st, r1st_analytic);
legend('Inverse crime', 'No crime', 'Analytic', 'location', 'southwest', 'FontSize',15)
xlabel('$\alpha$','Interpreter','Latex','FontSize',15)
ylabel('Relative error','FontSize',15)

%Find the corresponding alphas with the smallest relative errors.
[minr1_crime,idx1_crime]       = min(r1st_crime);     optalpha_1st_crime    = alpha_relerr_1st(idx1_crime);
[minr1_nocrime,idx1_nocrime]   = min(r1st_nocrime);   optalpha_1st_nocrime  = alpha_relerr_1st(idx1_nocrime);
[minr1_analytic,idx1_analytic] = min(r1st_analytic);  optalpha_1st_analytic = alpha_relerr_1st(idx1_analytic);

%Find the corresponding solutions
xopt1_crime    = TikhNN(A1,b1,optalpha_1st_crime, L);
xopt1_nocrime  = TikhNN(A2,b2,optalpha_1st_nocrime, L);
xopt1_analytic = TikhNN(A3,b3,optalpha_1st_analytic, L);

figure
semilogx(alpha_relerr_1st, r1st_crime);
hold on
semilogx(alpha_relerr_1st, r1st_nocrime);
hold on
semilogx(alpha_relerr_1st, r1st_analytic);
legend('Inverse crime', 'No crime', 'Analytic', 'location', 'southwest', 'FontSize',15)
xlabel('$\alpha$','Interpreter','Latex','FontSize',15)
ylabel('Relative error','FontSize',15)

figure
showDistribution(xopt1_crime,ginfo); title(sprintf('Inverse crime, relative error: %.4f', minr1_crime))

figure
showDistribution(xopt1_nocrime,ginfo); title(sprintf('No crime, relative error: %.4f', minr1_nocrime))

figure
showDistribution(xopt1_analytic,ginfo); title(sprintf('Analytic, relative error: %.4f', minr1_analytic))
