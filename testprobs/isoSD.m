function [A, b, x, L, ginfo] = isoSD(varargin)
%Constructs a isotropic slowing down test case.

% -- Default parameters for slowing down distribution --
%Parameters for the (vpara,vperp)-grid.
vparamin = -1.4e7;
vparamax = 1.4e7;
vperpmin = 1e5;
vperpmax = 1.4e7;
vparadim = 40; vperpdim = 20;

% -- Default parameters for A matrix --
phi = [10 20 40 70 85];
Emax = 4e6; udim = 100;

% -- Default right hand side
rhs = 'analytic';

% -- Unpack default parameters --
validvars = {'vparamin','vparamax','vperpmin','vperpmax', ...
             'vparadim','vperpdim','phi','umax','du'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end

% -- Construct the grid, x, u and A --
[vpara, vperp, ginfo] = construct_vgrid(vparadim,  vperpdim, ...
                                       'vperpmin', vperpmin, ...
                                       'vperpmax', vperpmax, ...
                                       'vparamin', vparamin, ...
                                       'vparamax', vparamax);

x = isoSDx(vpara,vperp);                                   
u = construct_uvec('Emax', 4e6, 'udim', 100);
A = transferMatrix(vpara,vperp,phi,u);

switch rhs
    case 'inverse crime'
        [b, e] = generate_noisy_b(A,xtrue); %Basically b = A*x + noise.
    case 'no crime'
        vparadimfine = 100; vperpdimfine = 50;
        [vparafine, vperpfine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
        xfine = isoSDx(vparafine,vperpfine);
        Afine = transferMatrix(vparafine,vperpfine,phi,u);
        [b, e] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.
    case 'analytic'
        gu = isoSDg(u,phi);
        du = u(2)-u(1);
        b = gu*du;
        [b, e] = add_noise(b);
    otherwise
        error('Unknown right hand side.')
end

%Normalize with this noise.
[A, b] = error_normalization(A,b,e); 

%Regularization matrix L
L = reguL(vpara, vperp);
end