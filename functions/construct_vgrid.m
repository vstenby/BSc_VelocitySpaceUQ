function [vpara, vperp, gridinfo] = construct_vgrid(vparamin, vparamax, vparadim, vperpmin, vperpmax, vperpdim)
%Construct VPERP and VPARA which are either vectors or matrices.

vpara_linspace = linspace(vparamin, vparamax, vparadim);
vperp_linspace = linspace(vperpmin, vperpmax, vperpdim);

%dvpara = (vparamax-vparamin)/(vparadim-1);
%dvperp = (vperpmax-vperpmin)/(vperpdim-1);

[vpara, vperp]=meshgrid(vpara_linspace,vperp_linspace);

gridinfo.vpara_ax = vpara_linspace;
gridinfo.vperp_ax = vperp_linspace;

end