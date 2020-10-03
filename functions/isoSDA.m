function W = isoSDW(ubroadening, isoSDFinfo, isoSDSinfo)
%Create the A matrix for the isotropic slowing down.

vpara = isoSDFinfo.vpara;
vperp = isoSDFinfo.vperp;
phi   = isoSDSinfo.phi;
u     = isoSDSinfo.u;
dures = isoSDSinfo.du * ubroadening;

W = transferMatrix(vpara, vperp, phi, u, dures);
end