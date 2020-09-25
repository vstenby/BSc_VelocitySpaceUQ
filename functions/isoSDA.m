function A = isoSDA(ubroadening, isoSDXinfo, isoSDSinfo)
%Create the A matrix for the isotropic slowing down.

vpara = isoSDXinfo.vpara;
vperp = isoSDXinfo.vperp;
phi   = isoSDSinfo.phi;
u     = isoSDSinfo.u;
dures = isoSDSinfo.du * ubroadening;

A = transferMatrix(vpara, vperp, phi, u, dures);
end