%% Drifting bi-Maxwellian
clear, clc, close all

% Setting the parameters.
Mp   = 1.6726e-27;
Qe   = 1.6021917e-19;
Mi   = 4*Mp;

%boundaries of the (E,p)-space
Emin = 10e3;
Emax = 4e6;

%plasma parameters, isotropic means Tpara=Tperp
Tpara=200e3; %eV
Tperp=200e3; %eV

ni=1e19;     %1/m^3
vparadrift=5e6;

%% Constructing the (vpara, vperp) grid.

vparamin = -1.4e7;
vparamax = 1.4e7;
vparadim = 100;
dvpara   = (vparamax-vparamin)/(vparadim-1);
vpara    = linspace(vparamin, vparamax, vparadim);

vperpmin = 1e5;
vperpmax = 1.4e7;
vperpdim = 50;
dvperp   = (vperpmax-vperpmin)/(vperpdim-1);
vperp    = linspace(vperpmin, vperpmax, vperpdim);

[VPARA,VPERP]=meshgrid(vpara,vperp);
 
%observation angles
phivec=[10 20 40 70 85];

%spectral data points
umin=sqrt(2*Emin*Qe/Mi);
umax=sqrt(2*Emax*Qe/Mi);

%number of points per spectrum
udim=200;

%spectral resolution of the measurements divided by bin width u of the spectra
ubroadening=3;

%% Using biMax to construct a drifting bi-Maxwellian

[X, S, A] = biMax(VPARA, VPERP, Tpara, Tperp, vparadrift, phivec, umin, umax, udim, ubroadening);

figure
sgtitle('2D bi-Maxwellian with drift')
subplot(1,2,1)
imagesc(X); colorbar(); axis xy; axis image; title('Density function')

subplot(1,2,2)
plot(S, '.')
hold on
plot(A*X(:)); 
legend('Analytic projections', 'Numeric projections')

%%

alpha = 100000;
x_an_clean = mosek_TikhNN(A, S, alpha);
x_an_noisy = mosek_TikhNN(A, add_noise(S), alpha);

x_nu_clean = mosek_TikhNN(A, A*X(:), alpha);
x_nu_noisy = mosek_TikhNN(A, add_noise(A*X(:)), alpha);

figure
imagesc(X); colorbar(); axis xy; axis image;

figure
imagesc(reshape(x_an_clean, vperpdim, vparadim)); colorbar(); axis xy; axis image;
title('TikhNN on clean, analytic S. alpha = 100')

figure
imagesc(reshape(x_an_noisy, vperpdim, vparadim)); colorbar(); axis xy; axis image;
title('TikhNN on noisy, analytic S. alpha = 100')

figure
imagesc(reshape(x_nu_clean, vperpdim, vparadim)); colorbar(); axis xy; axis image;
title('TikhNN on clean, numeric S. alpha = 100')

figure
imagesc(reshape(x_nu_noisy, vperpdim, vparadim)); colorbar(); axis xy; axis image;
title('TikhNN on noisy, numeric S. alpha = 100')

%% Normalizing first
%We normalize with our estimate of the error on the signal, 0.01. 
%This comes from the add noise function.

[S_an_norm, A_an_norm, factor_an] = normalize_problem(add_noise(S), A, 0.01); 
[S_nu_norm, A_nu_norm, factor_nu] = normalize_problem(add_noise(A*X(:)), A, 0.01);

alpha = 0.001;
x_an_noisy2 = mosek_TikhNN(A_an_norm, S_an_norm, alpha) * factor_an;
x_nu_noisy2 = mosek_TikhNN(A_nu_norm, S_nu_norm, alpha) * factor_nu;

figure
imagesc(reshape(x_an_noisy2, vperpdim, vparadim)); colorbar(); axis xy; axis image;
title('TikhNN on noisy, analytic S with normalizing before solving')

figure
imagesc(reshape(x_nu_noisy2, vperpdim, vparadim)); colorbar(); axis xy; axis image;
title('TikhNN on noisy, numeric S with normalizing before solving')







