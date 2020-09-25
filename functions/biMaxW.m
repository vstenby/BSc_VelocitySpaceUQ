function W = biMaxW(ubroadening, biMaxFinfo, biMaxSinfo)
%Create the W matrix for the bi-Maxwellian.

vpara = biMaxFinfo.vpara;
vperp = biMaxFinfo.vperp;
phi   = biMaxSinfo.phi;
u     = biMaxSinfo.u;
dures = biMaxSinfo.du * ubroadening;

W = transferMatrix(vpara, vperp, phi, u, dures);
end