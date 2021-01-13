function [A, b, x, L, ginfo] = biMax(varargin)
%Constructs a bi-Maxwellian testcase.

% -- Default parameters for bi-Maxellian --
%Grid parameters for true solution
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;
vparadim = 40; vperpdim = 20;

% -- Default parameters for A matrix --
phi = [10 20 40 70 85];
umax = 5e6; du = 1e5;

% -- Default right hand side --
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

x = biMaxx(vpara, vperp);
u = construct_uvec('umax',umax,'du',du);
A = transferMatrix(vpara,vperp,phi,u);

switch rhs
    case 'inverse crime'
        [b, e] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 
    case 'no crime'
        vparadimfine = 100; vperpdimfine = 50;
        [vparafine, vperpfine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
        xfine = biMaxx(vparafine, vperpfine); 
        Afine = transferMatrix(vparafine,vperpfine,phi,u); 
        [b, e] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.
    case 'analytic'
        gu = biMaxg(u,phi);
        du = u(2)-u(1);
        b = du*gu;
        [b, e] = add_noise(b);
    otherwise
        error('Unknown right hand side.');
end

%Normalize with this noise.
[A, b] = error_normalization(A,b,e);

%Regularisation matrix L
L = reguL(vpara,vperp);
end