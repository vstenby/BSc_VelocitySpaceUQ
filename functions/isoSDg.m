function [g, info] = isoSDg(u,phi,varargin)
% Evaluates the analytic projected velocity distribution of the
% slowing down fast-ion velocity distribution.
%
% Usage: 
%    ``g = isoSDg(u)``
%
%    ``g = isoSDg(u,phi)``
%
%    ``[g, info] = isoSDg(u,phi,varargin)``
%
% Inputs:
%    * **u**:               Vector of projected velocities.
%
%    * **phi**:             Observation angles.
%
% Optional inputs:
%    * **Ecrit**:               Critical energy. Default : ``44*20000 eV``
%   
%    * **Ebirth**:              Birth energy. Default : ``3.5e6 eV``
%   
%    * **Ebirthwidth**:         Birth energy width. Default : ``6e4``
%
%    * **Mi**:                  Mass of ions. Default : ``4*Mp``
%
%    * **ne**:                  Number of ions. Default : ``1e19``
%
% Output:
%    * **g**:               Projected velocity distribution.
%
%    * **info**:            MATLAB struct containing all inputs (including optional inputs) above.

%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

if nargin == 1
   %Since the isotropic slowing-down distribution doesn't depend on 
   %phi, we don't need it.
   phi = [];  
end

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
%Set default values for optional parameters.
Ecrit=44*20000;  %eV
Ebirth=3.5e6; %eV
Ebirthwidth=6e4;
Mi = 4*Mp;
ne = 1e19; 

%Unpack the varargin and evaluate.
validvars = {'Ecrit','Ebirth','Ebirthwidth','Mi','ne'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

%slowing-down
vcrit=sqrt(2*Ecrit*Qe/Mi);
vbirth=sqrt(2*Ebirth*Qe/Mi);

g=ne/((log(1+(vbirth/vcrit)^3))*4*vcrit)*(log(abs((vbirth^2-vcrit*vbirth+vcrit^2)./(u.^2-vcrit*abs(u)+vcrit^2)))+...
    log(((abs(u)+vcrit)/(vbirth+vcrit)).^2)+2*sqrt(3)*(atan((2*vbirth-vcrit)/(sqrt(3)*vcrit))-atan((2*abs(u)-vcrit)/(sqrt(3)*vcrit)))).*erfc((0.5*Mi*u.^2/Qe-Ebirth)/Ebirthwidth)/2; 

%Apply the heaviside function.
g = g.*heaviside(vbirth - abs(u));
g = reshape(g, length(g), 1);

if ~isempty(phi)
    g = repmat(g, [length(phi),1]); %Since b doesn't depend on phi, this is repeated.
end

%Saving the relevant parameters to the info structure.
info.u           = u;
info.phi         = phi;
info.Ecrit       = Ecrit;
info.Ebirth      = Ebirth;
info.Ebirthwidth = Ebirthwidth;
info.Mi          = Mi;
info.ne          = ne;
end

