function [g, info] = biMaxg(u,phi,varargin)
% Evaluates the analytic projected velocity distribution of the
% bi-Maxwellian velocity distribution.
%
% Usage: 
%    ``g = biMaxg(u)``
%
%    ``g = biMaxg(u,phi)``
%
%    ``[g, info] = biMaxg(u,phi,varargin)``
%
% Inputs:
%    * **u**:               Vector of projected velocities.
%
%    * **phi**:             Observation angles.
%
% Optional inputs:
%    * **Tpara**:           Parallel ion temperature. Default : ``20 keV``
%   
%    * **Tperp**:           Perpendicular ion temperature. Default : ``20 keV``
%   
%    * **vparadrift**:      Parallel-drift velocity. Default : ``5e5``
%
%    * **Mi**:              Mass of ion. Default : ``2*Mp``
%
%    * **ne**:              Number of ions. Default : ``1e19``
%
% Output:
%    * **g**:               Analytic projected velocity distribution.
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

g = [];
n_points = length(u);

%Constructing g analytically for each phi.
for i = 1:length(phi)
    phi_temp = phi(i);
    Teff=(Tperp*1000)*(sin(phi_temp/180*pi))^2+(Tpara*1000)*(cos(phi_temp/180*pi))^2; %Equation (72), denoted Tu.
    
    %We make length(u) points for each angle.
    idx1 = (i-1)*n_points + 1;
    idx2 = idx1+n_points - 1;
    
    %u_d is specified in Equation 73 but written inline here.
    g(idx1 : idx2) = ne*(Mi/(2*pi*Teff*Qe))^0.5*exp(-(Mi*(u-vparadrift*cos(phi_temp/180*pi)).^2)/(2*Teff*Qe)); %Equation (72)
end

%Return g as a column vector.
g = g'; 

%Saving the relevant parameters to the info structure.
info.u = u;
info.phi = phi;
info.Tpara = Tpara;
info.Tperp = Tperp;
info.vparadrift = vparadrift;
info.Mi = Mi;
info.ne = ne;
end