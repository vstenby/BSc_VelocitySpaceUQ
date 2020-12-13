%% Analysis of the relative error 
function biMax_sim_phicontour(sim_idx)
    %Add dependencies to path.
    addpath(genpath('../functions'))
    addpath(genpath('../../aux'))
    
    res = 0.25;
    phis = combvec(5:res:30, 30+res:res:55);
    phi3 = phis(1, sim_idx);
    phi4 = phis(2, sim_idx);
    
    r = min_relerr(phi3, phi4);
    
    save(sprintf('./biMax_sim_phicontour/r_%d.mat',sim_idx),'r')
    
end

function r = min_relerr(phi3, phi4)
phi  = [60 80];

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 35; vperpdim = 17;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

vparadimfine = 100; vperpdimfine = 50;
[vparafine, vperpfine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
xfine = biMaxx(vparafine, vperpfine); 

A = transferMatrix(vpara,vperp,[phi phi3 phi4],u);

Afine = transferMatrix(vparafine,vperpfine,[phi phi3 phi4],u); 
[b, e] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.
%Normalize with this noise.
[A, b] = error_normalization(A,b,e);   
%Regularization matrix L (1st order Tikhonov)
L = reguL(vperpdim, vparadim);
%A vector of alphas for finding the optimal regularization parameter.
alpha_relerr = logspace(-9,-4,100);
[~, ~, r1st] = TikhNN(A, b, alpha_relerr,  L, 'return_relerr', true, 'x_true', xtrue);
r = min(r1st);
end


