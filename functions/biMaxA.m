function A = biMaxA(ubroadening, biMaxxinfo, biMaxbinfo)
%Create the W matrix for the bi-Maxwellian.

vpara = biMaxxinfo.vpara;
vperp = biMaxxinfo.vperp;
phi   = biMaxbinfo.phi;
u     = biMaxbinfo.u;
dures = biMaxbinfo.du * ubroadening;

A = transferMatrix(vpara, vperp, phi, u, dures);
end