function [x, info] = isoSDx(vpara, vperp, varargin)
% Evaluates the isotropic slowing-down distribution on a (vpara,vperp)-grid.
%
% Usage: 
%    ``x = isoSDx(vpara,vperp)``
%
%    ``[x, info] = isoSDx(vpara,vperp,varargin)``
%
% Inputs:
%    * **vpara**:               Meshgrid for X to be evaluated on.
%
%    * **vperp**:               Meshgrid for X to be evaluated on.
%
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
%    * **x**:               The isotropic slowing-down distribution evaluated on the (vpara,vperp)-grid.
%
%    * **info**:            MATLAB struct containing all inputs (including optional inputs) above.


%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
%Set default values for optional inputs.
Ecrit = 44*20000;  %eV
Ebirth = 3.5e6;    %eV
Ebirthwidth = 6e4;
Mi  = 4*Mp;
ne  = 1e19;  
validvars = {'Ecrit','Ebirth','Ebirthwidth','Mi','ne'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

%slowing-down
vcrit=sqrt(2*Ecrit*Qe/Mi);
vbirth=sqrt(2*Ebirth*Qe/Mi);

fvpavpe3DSD = ne*3/(4*pi*log(1+(vbirth/vcrit)^3)).*heaviside(vbirth - sqrt(vpara.^2 + vperp.^2))./((vpara.^2+vperp.^2).^1.5 + vcrit^3);

x = fvpavpe3DSD.*(2*pi*vperp);

%Saving the relevant parameters to the info structure.
info.vpara       = vpara;
info.vperp       = vperp;
info.Ecrit       = Ecrit;
info.Ebirth      = Ebirth;
info.Ebirthwidth = Ebirthwidth;
info.Mi          = Mi;
info.ne          = ne;
end

