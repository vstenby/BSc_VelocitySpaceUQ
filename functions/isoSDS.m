function [S, info] = isoSDS(u, phi, Ecrit, Ebirth, Ebirthwidth, options)
% This function evaluates the isotropic slowing-down distribution.
%
%
% If not specified, these will be the default values:
% Tpara = 200e3 eV
% Tperp = 200e3 eV
% vparadrift = 5e6
% Mi = 4*(1.6726e-27)
% ne = 1e19
% nphi = 1


%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

%A bit of code duplication, perhaps ask Jakob how to do this better.
if nargin == 2
    Ecrit=44*20000;  %eV
    Ebirth=3.5e6; %eV
    Ebirthwidth=6e4;
    options.Mi = 4*Mp;
    options.ne = 1e19; 
elseif nargin == 4
    Ebirthwidth=6e4;
    options.Mi = 4*Mp;
    options.ne = 1e19; 
elseif nargin == 5
    %options not given.
    options.Mi = 4*Mp;
    options.ne = 1e19;
elseif nargin == 6
    %All options given.
    if ~isstruct(options), error('6th argument should be a struct'), end
    if ~isfield(options,'Mi'), options.Mi = 4*Mp; end %Default if not specified.
    if ~isfield(options,'ne'), options.ne = 1e19; end %Default if not specified.
else
    %Too many or too few inputs.
    error('Wrong number of inputs')
end

%Fetching other parameters from options.
ne = options.ne;
Mi = options.Mi;

%Make sure we get u as specified.
if isstruct(u)
    uvec = construct_uvec(u, Mi);
else
    %Assume u is a vector.
    uvec = u;
end

%slowing-down
vcrit=sqrt(2*Ecrit*Qe/Mi);
vbirth=sqrt(2*Ebirth*Qe/Mi);

guSD=ne/((log(1+(vbirth/vcrit)^3))*4*vcrit)*(log(abs((vbirth^2-vcrit*vbirth+vcrit^2)./(uvec.^2-vcrit*abs(uvec)+vcrit^2)))+...
    log(((abs(uvec)+vcrit)/(vbirth+vcrit)).^2)+2*sqrt(3)*(atan((2*vbirth-vcrit)/(sqrt(3)*vcrit))-atan((2*abs(uvec)-vcrit)/(sqrt(3)*vcrit)))).*erfc((0.5*Mi*uvec.^2/Qe-Ebirth)/Ebirthwidth)/2; 

guSD = guSD';

S = guSD;
S = repmat(S, [length(phi),1]); %Since S doesn't depend on phi, this is repeated.


%Saving the relevant parameters to the info structure.
info.u           = uvec;
info.du          = info.u(2)-info.u(1);
info.Ecrit       = Ecrit;
info.Ebirth      = Ebirth;
info.Ebirthwidth = Ebirthwidth;
info.Mi          = Mi;
info.ne          = ne;
info.phi         = phi;




end

