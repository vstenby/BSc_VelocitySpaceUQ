function [vpara, vperp, ginfo] = construct_vgrid(varargin)
%   Construct a (vpara,vperp)-grid. 
%
%   Usage:
%   [vpara, vperp, gridinfo] = construct_vgrid()
%   [vpara, vperp, gridinfo] = construct_vgrid(vparadim, vperpdim)
%   
%   Multiple arguments can be given in the following way:
%   construct_vgrid(vparadim,vperpdim, 'vperpmin', -1e8 ...)   
%   where you can pass 'vparamin', 'vparamax', 'vperpmin' and 'vperpmax'
%   as optional arguments.
%
%   Viktor Stenby Johansson, Fall 2020.

%Default parameters.
vparamin = -1.4e7;
vparamax = 1.4e7;
vparadim = 100;
vperpmin = 1e5;
vperpmax = 1.4e7;
vperpdim = 50;
    
switch nargin
    case 0
        vparadim = 100;         vperpdim = 50;
        default_parameters = 1;
    case 1
        vparadim = varargin{1}; vperpdim = 50;
        default_parameters = 1;
    case 2
        vparadim = varargin{1}; vperpdim = varargin{2};
    otherwise
        vparadim = varargin{1}; vperpdim = varargin{2};
        %Unpack the remaining varargin.
        optargin  = {varargin{3:end}};
        noptargin = length(optargin);
        if mod(noptargin,2) == 1 
            error('Wrong number of optional argin'); 
        else
            %Unpack the reamining argin.
            for i=1:2:noptargin
               arg = optargin{i};
               switch arg
                   case 'vparamin'
                       vparamin = optargin{i+1};
                   case 'vparamax'
                       vparamax = optargin{i+1};
                   case 'vperpmin'
                       vperpmin = optargin{i+1};
                   case 'vperpmax'
                       vperpmax = optargin{i+1};
                   otherwise
                       error('Unrecognized varargin');
               end
            end 
        end
end

vpara_linspace = linspace(vparamin, vparamax, vparadim);
vperp_linspace = linspace(vperpmin, vperpmax, vperpdim);

[vpara, vperp] = meshgrid(vpara_linspace,vperp_linspace);

ginfo.vparamin = vparamin;
ginfo.vparamax = vparamax;
ginfo.vparadim = vparadim;
ginfo.vpara_ax = vpara_linspace;

ginfo.vperpmin = vperpmin;
ginfo.vperpmax = vperpmax;
ginfo.vperpdim = vperpdim;
ginfo.vperp_ax = vperp_linspace;
end