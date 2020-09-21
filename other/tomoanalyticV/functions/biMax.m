function [X, S, A] = biMax(vpara, vperp, Tpara, Tperp, vparadrift, phivec, ...
                           umin, umax, udim, ubroadening)         
%   This function calculates the analytic density function X, 
%   the analytic projection S, and the forward model matrix A, such that
%   A*X(:) is the numeric projection S.
%   
%   s174483@student.dtu.dk
%
%
%   Input:
%   
%   vperp
%   vpara
%   Tpara
%   Tperp
%   vparadrift
%   phivec
%
%
%   Output:
%   
%   X: Density function for the parallel-drifting bi-Maxwellian distribution.
%   S: Projection for the parallel-drifting bi-Maxwellian distribution.


%   Relevant papers:
%
%   [A]
%
%   Bi-Maxwellian, slowing-down, and ring
%   velocity distributions of fast ions in
%   magnetized plasmas

%TODO: Should vperpdrift be included?
%      Should 3D be implemented as well?

Mp   = 1.6726e-27;
Qe   = 1.6021917e-19;
Mi   = 4*Mp;
ne = 1e19;

%define thermal velocities
vthpara=sqrt(2*Tpara*Qe/Mi);
vthperp=sqrt(2*Tperp*Qe/Mi);

%Equation 69 in [A].
X = (2*ne*vperp)/(sqrt(pi)*vthpara*vthperp.^2) .* exp(-((vpara-vparadrift)/vthpara).^2-(vperp/vthperp).^2);

%Check that vperp and vpara (often grids) have the same size.
[x1, x2] = size(vperp);
[y1, y2] = size(vpara);

assert(x1 == y1 & x2 == y2, 'Grid sizes do not match.');


%Used for constructing A as well.
du=(umax-umin)/(udim/2-1);
dures=ubroadening*du;
uvec = linspace(-umax, umax, udim); %Why do we use -umax here?


S = [];
n_points = length(uvec);

%Constructing S analytically.
for i = 1:length(phivec)
    phi = phivec(i);
    Teff=Tperp*(sin(phi/180*pi))^2+Tpara*(cos(phi/180*pi))^2; %Equation (72), denoted Tu.
    
    %We make length(uvec) points for each angle.
    idx1 = (i-1)*n_points + 1;
    idx2 = idx1+n_points - 1;
    
    %u_d is specified in Equation 73 but written inline here.
    S(idx1 : idx2) = ne*(Mi/(2*pi*Teff*Qe))^0.5*exp(-(Mi*(uvec-vparadrift*cos(phi/180*pi)).^2)/(2*Teff*Qe)); %Equation (72)
end

%Return the projection S as a column-vector.
S = S';

A = transfer_matrix(vpara, vperp, phivec, uvec, dures);
end