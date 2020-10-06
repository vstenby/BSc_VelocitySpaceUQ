%% Investigate the number of samples needed.
%This script is investigating the number of samples.
clear, clc, close all

%Set up the default problem.
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

%Observation angles
phi=[10 20 40 70 85];

[b, binfo] = biMaxb(ustruct,phi);

A = biMaxA(ubroadening, xinfo, binfo);






