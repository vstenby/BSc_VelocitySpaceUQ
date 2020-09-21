function [S, info] = biMaxS(u, phi, Tpara, Tperp, vparadrift, options)
% This function returns the analytical projections of the bi-Maxwellian
% distribution function.
%
% 
% If not specified, these will be the default values:
% Tpara = 200e3 eV
% Tperp = 200e3 eV
% vparadrift = 5e6
% Mi = 4*(1.6726e-27)
% ne = 1e19

%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

%A bit of code duplication, perhaps ask Jakob how to do this better.
if nargin == 2
    %uvec phivec specified,
    Tpara=200e3; %eV
    Tperp=200e3; %eV
    vparadrift = 5e6; 
    options.Mi = 4*Mp;
    options.ne = 1e19;  
elseif nargin == 4
    %tpara tperp specified
    vparadrift = 5e6; 
    options.Mi = 4*Mp;
    options.ne = 1e19;  
elseif nargin == 5
    %vparadrift specified
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

S = [];
n_points = length(uvec);
phivec = phi;

%Constructing S analytically.
for i = 1:length(phivec)
    phi = phivec(i);
    Teff=Tperp*(sin(phi/180*pi))^2+Tpara*(cos(phi/180*pi))^2; %Equation (72), denoted Tu.
    
    %We make length(uvec) points for each angle.
    idx1 = (i-1)*n_points + 1;
    idx2 = idx1+n_points - 1;
    
    %u_d is specified in Equation 73 but written inline here.
    S(idx1 : idx2) = ne*(Mi/(2*pi*Teff*Qe))^0.5*exp(-(Mi*(uvec-vparadrift*cos(phi/180*pi)).^2)/(2*Teff*Qe)); %Equation (72)
end

%Return S as a column vector.
S = S'; 

%Saving the relevant parameters to the info structure.
info.u = uvec;
info.du = info.u(2)-info.u(1);
info.phi = phivec;
info.Tpara = Tpara;
info.Tperp = Tperp;
info.vparadrift = vparadrift;
info.Mi = Mi;
info.ne = ne;

end