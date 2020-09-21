function [X, S, A] = isoSD(vpara, vperp, phivec, umin, umax, udim, ubroadening)

Mp   = 1.6726e-27;
Qe   = 1.6021917e-19;
Mi   = 4*Mp;
ne = 1e19; 

%slowing-down
Ecrit=44*20000;  %eV
vcrit=sqrt(2*Ecrit*Qe/Mi);
Ebirth=3.5e6; %eV
vbirth=sqrt(2*Ebirth*Qe/Mi);
Ebirthwidth=6e4;

%Used for constructing A as well.
du=(umax-umin)/(udim/2-1);
dures=ubroadening*du;
uvec = linspace(-umax, umax, udim); %Why do we use -umax here?

fvpavpe3DSD=ne*3/(4*pi)/(log(1+(vbirth/vcrit)^3))./((vpara.^2+vperp.^2).^1.5+vcrit^3).*erfc((0.5*Mi*(vpara.^2+vperp.^2)/Qe-Ebirth)/Ebirthwidth)/2;
X = fvpavpe3DSD.*(2*pi*vperp);

guSD=ne/((log(1+(vbirth/vcrit)^3))*4*vcrit)*(log(abs((vbirth^2-vcrit*vbirth+vcrit^2)./(uvec.^2-vcrit*abs(uvec)+vcrit^2)))+...
    log(((abs(uvec)+vcrit)/(vbirth+vcrit)).^2)+2*sqrt(3)*(atan((2*vbirth-vcrit)/(sqrt(3)*vcrit))-atan((2*abs(uvec)-vcrit)/(sqrt(3)*vcrit)))).*erfc((0.5*Mi*uvec.^2/Qe-Ebirth)/Ebirthwidth)/2; 

guSD = guSD';

S = guSD;
S = repmat(S, [length(phivec),1]); %Since S doesn't depend on phi, this is repeated.

A = transfer_matrix(vpara, vperp, phivec, uvec, dures);


end