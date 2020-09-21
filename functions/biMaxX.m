function X = biMaxX(vpara, vperp, Tpara, Tperp, vparadrift, options)
%   DESCRIPTION GOES HERE!
%
%
%
%
%

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
X = (2*ne*vperp)/(sqrt(pi)*vthpara*vthperp.^2) .* exp(-((vpara-vparadrift)/vthpara).^2-(vperp/vthperp).^2);

end

