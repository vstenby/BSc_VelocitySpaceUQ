function [F, info] = biMaxF(vpara, vperp, Tpara, Tperp, vparadrift, options)
% This function evaluates the drifting bi-Maxwellian.
%
%
% If not specified, these will be the default values:
% Tpara = 200e3 eV - perhaps rewrite to 2e5 eV.?
% Tperp = 200e3 eV - perhaps rewrite to 2e5 eV.?
% vparadrift = 5e6
% Mi = 4*(1.6726e-27)
% ne = 1e19


%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

%A bit of code duplication, perhaps ask Jakob how to do this better.
if nargin == 2
    %Tpara, Tperp, vparadrift and options not given.
    Tpara=200e3; %eV
    Tperp=200e3; %eV
    vparadrift = 5e6; 
    options.Mi = 4*Mp;
    options.ne = 1e19;   
elseif nargin == 4
    %vparadrift and options not given.
    vparadrift =5e6; 
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

%Thermal velocities
vthpara=sqrt(2*Tpara*Qe/Mi);
vthperp=sqrt(2*Tperp*Qe/Mi);

%Equation 69.
F = (2*ne*vperp)/(sqrt(pi)*vthpara*vthperp.^2) .* exp(-((vpara-vparadrift)/vthpara).^2-(vperp/vthperp).^2);

if any(size(F) == 1)
    F = reshape(F, numel(F), 1);
end

%Saving the relevant parameters to the info structure.
info.vpara = vpara;
info.vperp = vperp;
info.Tpara = Tpara;
info.Tperp = Tperp;
info.vparadrift = vparadrift;
info.Mi = Mi;
info.ne = ne;
end

