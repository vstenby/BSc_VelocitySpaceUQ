%plot a bi-Maxwellian distribution ignoring potential energy
clear; clc

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

%slowing-down
Ecrit=44*20000;  %eV
vcrit=sqrt(2*Ecrit*Qe/Mi);
Ebirth=3.5e6; %eV
vbirth=sqrt(2*Ebirth*Qe/Mi);
Ebirthwidth=6e4;

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

%% 2D bi-Maxwellian with drift.

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

x = mosek_TikhNN(A, A*X(:), 0.1);

figure
imagesc(reshape(x, vperpdim, vparadim)); colorbar()

%% 2D Isotropic Slowing Down

[X, S, A] = isoSD(VPARA, VPERP, phivec, umin, umax, udim, ubroadening);

x = mosek_TikhNN(A, A*X(:), 0.1);

figure
subplot(1,2,1)
imagesc(reshape(X, vperpdim, vparadim)); colorbar()
subplot(1,2,2)
imagesc(reshape(x, vperpdim, vparadim)); colorbar()

