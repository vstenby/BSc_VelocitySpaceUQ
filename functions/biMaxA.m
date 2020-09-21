function A = biMaxA(ubroadening, biMaxXinfo, biMaxSinfo)
%Create the A matrix for the bi-Maxwellian.

vpara = biMaxXinfo.vpara;
vperp = biMaxXinfo.vperp;
phi   = biMaxSinfo.phi;
u     = biMaxSinfo.u;
dures = biMaxSinfo.du * ubroadening;

A = transferMatrix(vpara, vperp, phi, u, dures);
end