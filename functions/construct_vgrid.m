function [vpara, vperp, gridinfo] = construct_vgrid(vparamin, vparamax, vparadim, vperpmin, vperpmax, vperpdim)

if nargin == 0
    %Default values for construct_vgrid();
    vparamin = -1.4e7;
    vparamax = 1.4e7;
    vparadim = 100;
    vperpmin = 1e5;
    vperpmax = 1.4e7;
    vperpdim = 50;
end

vpara_linspace = linspace(vparamin, vparamax, vparadim);
vperp_linspace = linspace(vperpmin, vperpmax, vperpdim);

%dvpara = (vparamax-vparamin)/(vparadim-1);
%dvperp = (vperpmax-vperpmin)/(vperpdim-1);

[vpara, vperp]=meshgrid(vpara_linspace,vperp_linspace);

gridinfo.vpara_ax = vpara_linspace;
gridinfo.vperp_ax = vperp_linspace;

end