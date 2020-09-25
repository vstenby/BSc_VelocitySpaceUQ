function [X, info] = isoSDX(vpara, vperp, Ecrit, Ebirth, Ebirthwidth, options)
% This function evaluates the isotropic slowing-down distribution.
%
%
% If not specified, these will be the default values:
% Ecrit=44*20000;  %eV
% Ebirth=3.5e6; %eV
% Ebirthwidth=6e4;
% Mi = 4*(1.6726e-27)
% ne = 1e19

%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

%A bit of code duplication, perhaps ask Jakob how to do this better.
if nargin == 2
    Ecrit       =44*20000;  %eV
    Ebirth      = 3.5e6;    %eV
    Ebirthwidth = 6e4;
    options.Mi  = 4*Mp;
    options.ne  = 1e19;   
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

%slowing-down
vcrit=sqrt(2*Ecrit*Qe/Mi);
vbirth=sqrt(2*Ebirth*Qe/Mi);

fvpavpe3DSD=ne*3/(4*pi)/(log(1+(vbirth/vcrit)^3))./((vpara.^2+vperp.^2).^1.5+vcrit^3).*erfc((0.5*Mi*(vpara.^2+vperp.^2)/Qe-Ebirth)/Ebirthwidth)/2;
X = fvpavpe3DSD.*(2*pi*vperp);

%Saving the relevant parameters to the info structure.
info.vpara       = vpara;
info.vperp       = vperp;
info.Ecrit       = Ecrit;
info.Ebirth      = Ebirth;
info.Ebirthwidth = Ebirthwidth;
info.Mi          = Mi;
info.ne          = ne;

end

