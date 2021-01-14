function [x, info] = biMaxx(vpara, vperp, varargin)
% Evaluates the bi-Maxwellian fast-ion velocity
% distribution on a (vpara,vperp)-grid.
%
% Usage: 
%    ``x = biMaxx(vpara,vperp)``
%
%    ``[x, info] = biMaxx(vpara,vperp,varargin)``
%
% Inputs:
%    * **vpara**:               Meshgrid for X to be evaluated on.
%
%    * **vperp**:               Meshgrid for X to be evaluated on.
%
%
% Optional inputs:
%    * **Tpara**:               Parallel ion temperature.      Default : ``20 keV``
%   
%    * **Tperp**:               Perpendicular ion temperature. Default : ``20 keV``
%   
%    * **vparadrift**:          Parallel-drift velovity. Default : ``5e5 m/s``
%
%    * **Mi**:                  Mass of ion. Default : ``2*Mp``
%
%    * **ne**:                  Number of ions. Default : ``1e19``
%
% Output:
%    * **x**:               The Tikhonov solution for a given alpha.
%
%    * **info**:            MATLAB struct containing all inputs (including optional inputs) above.

%Physical constants
Mp = 1.6726e-27;    %Mass of proton
Qe = 1.6021917e-19; %Elementary charge

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
%Set default values for optional inputs.
Tpara = 20; %keV
Tperp = 20; %keV
vparadrift = 5e5;
Mi = 2*Mp;
ne = 1e19;

validvars = {'Tpara','Tperp','vparadrift','Mi','ne'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

vthpara=sqrt(2*Tpara*Qe*1000/Mi);
vthperp=sqrt(2*Tperp*Qe*1000/Mi);

%Equation 69.
x = (2*ne*vperp)/(sqrt(pi)*vthpara*vthperp.^2) .* exp(-((vpara-vparadrift)/vthpara).^2-(vperp/vthperp).^2);

if any(size(x) == 1)
    x = reshape(x, numel(x), x);
end

%Saving the relevant parameters to the info structure.
info.vpara = vpara;
info.vperp = vperp;
info.Tpara = Tpara;
info.Tperp = Tperp;
info.vparadrift = vparadrift;
info.Mi = Mi;
info.ne = ne;
end

