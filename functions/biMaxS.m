function S = biMaxS(u, phi, Tpara, Tperp, vparadrift, options)
%   DESCRIPTION GOES HERE!
% 
% u   : Either vector of values or structure
% phi : Vector 
%
%

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

%Make sure we get u as specified.
if isstruct(u)
    uvec = construct_uvec(ustruct);
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
end

%Auxil function to construct uvec.
function uvec = construct_uvec(ustruct)
umin = []; umax = []; du = []; udim = [];

if isfield(ustruct,'umin'), umin = ustruct.umin; end
if isfield(ustruct,'umax'), umax = ustruct.umax; end
if isfield(ustruct,'udim'), udim = ustruct.udim; end
if isfield(ustruct,'du'), du = ustruct.du; end

if ~isempty(udim) && ~isempty(du)
    error('udim and du cannot be specified simultaneously.')
elseif isempty(umin) && ~isempty(du)
    warning('du specified but not used.')
elseif isempty(umin) && isempty(udim)
    error('If umin is not specified, then udim should be.')
elseif isempty(umax)
    error('umax not specified')
end

if isempty(umin)
    %uvec should go from -umax to umax with udim points.
    uvec = linspace(-umax, umax, udim); 
elseif ~isempty(umin) && ~isempty(du)
    %uvec should look like this
    uvec = [-umax : du : -umin umin : du : umax];
end

end