function S = biMaxS(uvec, phivec, Tpara, Tperp, vparadrift, options)
%   DESCRIPTION GOES HERE!
%
%
%
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





end