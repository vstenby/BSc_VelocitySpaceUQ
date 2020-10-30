clear, clc, close all

crime = 1;

%set the grids. dvfine for computation of spectra, dv for inversion
%(vpara,vperp) grid
vparamax=1.4e7;
vparamin=-1.4e7;
vperpmin=1e5;
vperpmax=1.4e7;

%Observation angles
phi=[10 20 40 70 85];
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the fine forward model.
vparadimfine = 100; vperpdimfine = 50;
[vparafine, vperpfine, ginfofine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
xfine = isoSDx(vparafine, vperpfine); 

%Construct our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = isoSDx(vpara,vperp);

L = reguL(vperpdim, vparadim);

A = transferMatrix(vpara,vperp,phi,u); 
Afine = transferMatrix(vparafine,vperpfine,phi,u); 

if crime == 1
    [b, e] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 
    [A, b] = error_normalization(A,b,e);    %Normalize with this noise.
else
    [b, e] = generate_noisy_b(Afine,xfine); %Basically b = A*x + noise 
    [A, b] = error_normalization(A,b,e);    %Normalize with this noise.
end

alpha_relerr = logspace(-4, -2, 100);
%Find the relative error for 0th order and 1st order.
[~, ~, r0th] = TikhNN(A, b, alpha_relerr, [], 'return_relerr', true, 'x_true', xtrue);
[~, ~, r1st] = TikhNN(A, b, alpha_relerr,  L, 'return_relerr', true, 'x_true', xtrue);

semilogx(alpha_relerr, r0th)
hold on
semilogx(alpha_relerr, r1st)

%Calculate the optimal reconstructions for 
[minr0,idx0] = min(r0th); optalpha_0th = alpha_relerr(idx0);
[minr1,idx1] = min(r1st); optalpha_1st = alpha_relerr(idx1);
