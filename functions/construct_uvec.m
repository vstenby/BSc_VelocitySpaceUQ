function u = construct_uvec(varargin)
%   Viktor Stenby, s174483

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

for i=1:2:nargin
   arg = varargin{i}; 
   switch arg
       case 'umin'
           umin = varargin{i+1};
       case 'umax'
           umax = varargin{i+1};
       case 'Emin'
           Emin = varargin{i+1};   
       case 'Emax'
           Emax = varargin{i+1};
       case 'udim'
           udim = varargin{i+1};
       case 'du'
           du   = varargin{i+1};
       otherwise
           error('Wrong argument.');
   end      
end

if ~exist('Mi','var')
    %Default value of Mi.
    Mi = 4*Mp;
end

% If Emin and Emax are given, then these should be converted to 
% umin and umax. 

if exist('Emin','var')
    umin=sqrt(2*Emin*Qe/Mi);
end

if exist('Emax','var')
    umax=sqrt(2*Emax*Qe/Mi);
end

% By now, umin and umax should either be specified or set to default.
if ~exist('umax','var')
   %umax is set to default
   umax = sqrt(2*4e6*Qe/Mi);
end

if ~exist('umin','var')
   %umin is set to default 
   umin = sqrt(2*10e3*Qe/Mi);
end


if isempty(umin) %if umin = [] by the user.
   if ~exist('udim','var') && ~exist('du','var')
       %umin = [], neither udim or du are specified
       u = linspace(-umax,umax,200); 
   elseif exist('du', 'var')
       %umin = [], du specified
       u = -umax:du:umax;
   elseif exist('udim', 'var')
       %umin = [], udim specified
       u = linspace(-umax,umax,udim);
   else
       error('Something went wrong.')
   end
end

if ~isempty(umin) && exist('du','var')
    %u should look like this.
    u = [-umax : du : -umin umin : du : umax];
end

if ~isempty(umin)
   if exist('udim','var')
       %udim should be an even number.
       if mod(udim,2) == 1
           error('udim should be an even number')
       else
           u = [linspace(-umax, -umin, udim/2), linspace(umin, umax, udim/2)];
       end
   else
      u = [linspace(-umax, -umin, 100), linspace(umin, umax, 100)]; 
   end
end
end


function result = isGiven(arg,args)
result = any(strcmpi(arg,args));
end