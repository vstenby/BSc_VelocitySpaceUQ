function [vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim, varargin)
% Function for constructing (vpara,vperp)-grid on which we can evaluate our
% distribution functions.
%
% Usage: 
%    ``[vpara, vperp] = construct_vgrid(vparadim, vperpdim)``
%
%    ``[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim, varargin)`` 
%
% Inputs:
%    * **vparadim**:        The dimensions in the parallel dimension.
%
%    * **vperpdim**:        The dimensions in the perpendicular dimension.
%
% Optional inputs:
%    * **vparamin**:        Write description here.
%
%    * **vparamax**:        Write description here. 
%
%    * **vparadim**:        Write description here.
%
%    * **vperpmin**:        Write description here
%
%    * **vperpmax**:        Write description here.
%
%    * **vperpdim**:        Write description here.
%
% Output:
%    * **vpara**:           Write description here.
%
%    * **vperp**:            Write description here.
%
%    * **ginfo**:            Write description here.


switch nargin
    case 0
        vparadim = 100; vperpdim = 50;
    case 1
        vperpdim = 50;
end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
% Default values of optional inputs
vparamin = -1.4e7;
vparamax = 1.4e7;

vperpmin = 1e5;
vperpmax = 1.4e7;

validvars = {'vparamin','vparamax','vperpmin','vperpmax'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - -

if ~exist('vparadim','var') || ~exist('vperpdim','var')
    error('vparadim or vperpdim is not defined.')
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