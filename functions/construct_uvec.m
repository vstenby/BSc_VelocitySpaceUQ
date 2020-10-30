function u = construct_uvec(varargin)
% Construct the u vector.
%
% Usage: 
%    ``u = construct_uvec()`` 
%
%    ``u = construct_uvec(varargin)`` 
%
% Optional inputs:
%    * **Emin**:        Write description here.
%
%    * **Emax**:        Write description here.
%
%    * **umin**:        Write description here. 
%
%    * **umin**:        Write description here. 
%
%    * **du**:          Write description here.
%
%    * **udim**:        Write description here.
%
% Output:
%    * **u**:           Write some description here.

Qe = 1.6021917e-19; %Elementary charge
Mp = 1.6726e-27;    %Mass of proton

%Unpack varargin
if mod(nargin,2) == 1; error('Odd number of inputs'); end
args = {varargin{1:2:end}};

% -- Checks for proper usage --
if isGiven('umin',args) && isGiven('Emin',args)
    error('Both umin and Emin specified.')
end

if isGiven('umax',args) && isGiven('Emax',args)
    error('Both umax and Emax specified.')
end

if isGiven('du',args) && isGiven('udim',args)
   error('Both du and udim specified.') 
end
% -- End of checks --

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
% Default values of Mi
Mi   = 4*Mp;
udim = 200; %Default udim.

%Read varargin.
validvars = {'umin','umax','Emin','Emax','udim','du','Mi'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end

%If du is given as an optional input, we should clear udim.
if exist('du','var'), clear udim, end

if exist('Emin','var') && ~exist('umin','var')
%We are given Emin and not umin, convert Emin to umin.
umin=sqrt(2*Emin*Qe/Mi);
elseif ~exist('Emin','var') && exist('umin','var')
%We are not given Emin, should be fine as it is.
elseif exist('Emin','var') && exist('umin','var')
%We are given both, throw an error.
error('Both Emin and umin are given')
elseif ~exist('Emin','var') && ~exist('umin','var')
%We are given neither, which is fine.
end

if exist('Emax','var') && ~exist('umax','var')
%We are only given Emax, calculate umax.
umax = sqrt(2*Emax*Qe/Mi);
elseif ~exist('Emax','var') && exist('umax','var')
%We are not given Emin, should be fine as it is.
elseif exist('Emax','var') && exist('umax','var')
%We are given both, throw an error.
error('Both Emax and umax are given')
elseif ~exist('Emax','var') && ~exist('umax','var')
%We are given neither, use default umax.
umax = sqrt(2*4e6*Qe/Mi); %default Emax of 4e6
end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

% If Emin and Emax are given, then these should be converted to 
% umin and umax. 

if ~exist('umin','var')
   if exist('udim','var')
       u = linspace(-umax, umax, udim);
   elseif exist('du','var')
       u = -umax:du:umax;
   else
       error('Missing udim/du')
   end
end

if exist('umin','var')
   if exist('udim','var')
       if mod(udim,2) == 1, error('udim should be an even number.'), end
       u = [linspace(-umax, -umin, udim/2), linspace(umin, umax, udim/2)];
   elseif exist('du','var')
       u = [-umax : du : -umin umin : du : umax];
   else
       error('Missing udim/du')
   end
end
end


function result = isGiven(arg,args)
result = any(strcmpi(arg,args));
end