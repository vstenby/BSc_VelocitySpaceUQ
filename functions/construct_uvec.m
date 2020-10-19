function uvec = construct_uvec(ustruct, Mi)
%Auxil function to construct uvec.

Qe = 1.6021917e-19; %Elementary charge
if nargin < 2
    Mp = 1.6726e-27;    %Mass of proton
    Mi = 4*Mp;
end

umin = []; umax = []; du = []; udim = [];


if isfield(ustruct,'umin'), umin = ustruct.umin; end
if isfield(ustruct,'umax'), umax = ustruct.umax; end
if isfield(ustruct,'udim'), udim = ustruct.udim; end
if isfield(ustruct,'du'), du = ustruct.du; end
if isfield(ustruct,'Emin')
    if ~isempty(umin), warning('umin and Emin specified, umin calculated using Emin.'), end
    umin=sqrt(2*ustruct.Emin*Qe/Mi);
end
if isfield(ustruct,'Emax')
    if ~isempty(umax), warning('umax and Emax specified, umax calculated using Emax.'), end
    umax=sqrt(2*ustruct.Emax*Qe/Mi);
end

if ~isempty(udim) && ~isempty(du)
    error('udim and du cannot be specified simultaneously.')
elseif ~isempty(udim) && ~isempty(du)
    warning('du specified but not used.')
elseif isempty(umin) && (isempty(udim) && isempty(du))
    error('If umin is not specified, then udim or du should be.')
elseif isempty(umax)
    error('umax/Emax is not specified.')
elseif isempty(udim) && isempty(du)
    error('Neither udim or du are specified.')
end

if isempty(umin)
    if ~isempty(udim)
        %uvec should go from -umax to umax with udim points.
        uvec = linspace(-umax, umax, udim); 
    elseif ~isempty(du)
        uvec = -umax:du:umax;
    else
        error('Something went wrong.')
    end
elseif ~isempty(umin) && ~isempty(du)
    %uvec should look like this
    uvec = [-umax : du : -umin umin : du : umax];
elseif ~isempty(umin) && ~isempty(udim)
    %uvec should look like this
    if mod(udim,2)
        error('If umin is specified, then udim should be an even number.')
    end
    uvec = [linspace(-umax, -umin, udim/2), linspace(umin, umax, udim/2)];
else
    error('Something went wrong')
end

end