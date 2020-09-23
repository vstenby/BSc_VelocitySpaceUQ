%% Drifting bi-Maxwellian UQ
clear, clc, close all

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate the bi-Maxwellian on this grid with default values.
[X, biMaxXinfo] = biMaxX(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, biMaxSinfo] = biMaxS(ustruct,phi);
[S_noisy, s] = add_noise(S,0.01);

%Generate A from vpara, vperp, u and phi from biMaxX and biMaxS.
%First argument is ubroadening, which is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
A = biMaxA(3,biMaxXinfo,biMaxSinfo);

%Normalize by the measurement uncertainty
[Shat, Ahat] = measurement_normalization(S,A,s);

%Normalize with the 2-norm
[S2hat, A2hat, factor] = numeric_normalization(Shat, Ahat);

%% Setting up the sampler.
alpha = 3.4e-16;
xalpha = mosek_TikhNN(A2hat,S2hat,alpha);
[M, N] = size(A2hat);

%Initialize lambda0 and delta0.

lambda_est          = 1/norm(S2hat(:)-A2hat*xalpha(:))^2;
delta_est           = alpha*lambda_est;

% Initialization before sampling.  
nsamps     = 100;
lamsamp    = zeros(nsamps,1); lamsamp(1) = lambda_est;
delsamp    = zeros(nsamps,1); delsamp(1) = delta_est;
xsamp      = zeros(N,nsamps);
xtemp      = xalpha; 
xsamp(:,1) = xtemp(:);

%Hyperpriors
a0         = 1; 
t0         = 0.0001; 
a1         = 1; 
t1         = 0.0001;

for i=1:nsamps-1
   disp(i)
   Axtemp       = A2hat*xtemp; 
   lamsamp(i+1) = gamrnd(a0+M/2,1./(t0+norm(Axtemp(:)-S2hat(:))^2/2));
   Lxtemp       = xtemp; %L is the identity.
   delsamp(i+1) = gamrnd(a1+N/2,1./(t1+xtemp(:)'*Lxtemp(:)/2));
   alpha_temp   = delsamp(i+1)/lamsamp(i+1);
   x = mosek_TikhNN(A2hat, S2hat, alpha_temp);
end

%% Plot the chains

figure(1)
subplot(1,2,1)
plot(lamsamp)
subplot(1,2,2)
plot(delsamp)

figure(2)
showDistribution(x,gridinfo)
