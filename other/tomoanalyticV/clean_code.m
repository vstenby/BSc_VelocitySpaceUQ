% The relevant articles are:   
%
%   Bi-Maxwellian, slowing-down, and ring
%   velocity distributions of fast ions in
%   magnetized plasmas [A]
%
% and
%
%   On velocity-space sensitivity of fast-ion
%   D-alpha spectroscopy [B]
%
% which is available on DTU Findit.

%% Part 1 - Generating X.

%plot a bi-Maxwellian distribution ignoring potential energy
clear;

%No idea where these parameters come from.
Mp   = 1.6726e-27; 
Qe   = 1.6021917e-19;
MD=2*Mp;
Mi=4*Mp;


%Plasma parameters. Isotropic means Tpara=Tperp
%In [A], these are said to be ITER-relevant parameters.

Tpara=20;   %keV
Tperp=20;   %keV

%Very different values for Tpara and Tperp for the two examples.
%Tpara=200e3; %eV
%Tperp=200e3; %eV


ne=1e19;    %1/m^3  <--- this is said to be 1e20, not 1e19. Not sure if it matters.
ni=1e19;    %1/m^3
vdrift=5e5; %m/s  <--- this is mentioned in Figure 2 description box [A].
vparadrift=5e6; %m/s

%Thermal velocities
vthpara=sqrt(2*Tpara*Qe*1000/MD);
vthperp=sqrt(2*Tperp*Qe*1000/MD);

%set the grids. dv for computation of spectra, dv for inversion
%(vpara,vperp) grid
vparamax=4e6;
vparamin=-4e6;
vperpmin=1e4;
vperpmax=4e6;

%(vpara,vperp) grid
vparadim=40;
vperpdim=20;

%If we use the finer grid, then this takes ages.
%vparadim=100;
%vperpdim=50;

dvpara=(vparamax-vparamin)/(vparadim-1);
dvperp=(vperpmax-vperpmin)/(vperpdim-1);

vpara=vparamin:dvpara:vparamax;
vperp=vperpmin:dvperp:vperpmax;
[VPARA,VPERP]=meshgrid(vpara,vperp);
[rows,columns]=size(VPERP);

%bi-Maxwellian distribution function as slice in 3D
fvpavpe3DbiMax=ne/(pi^(3/2)*vthpara*vthperp^2)*exp(-((VPARA-vdrift)/vthpara).^2-(VPERP/vthperp).^2);

%fvpavpe in 2D
fvpavpe2DbiMax=fvpavpe3DbiMax.*(2*pi*VPERP);

%observation angles
phivecr=[10 20 40 70 85];

%This is equivalent to figure 4 in Bi-Maxwellian document.
imagesc(vpara/1e6,vperp/1e6,reshape(fvpavpe2DbiMax,vperpdim,length(vpara)))

%This is our X!
X = fvpavpe2DbiMax;

%Drifting bi-Maxwellian.
%My own attempt to implement Equation 69 in [A]. 
x_test = (2*ne*VPERP)/(sqrt(pi)*vthpara*vthperp.^2) .* exp(-((VPARA-vdrift)/vthpara).^2-(VPERP/vthperp).^2);

%Plotting the parallel-drifting bi-Maxwellian
f = @(x,y) (2*ne*y)/(sqrt(pi)*vthpara*vthperp.^2) .* exp(-((x-vdrift)/vthpara).^2-(y/vthperp).^2);
fsurf(f, [vparamin vparamax vperpmin vperpmax])

%% Part 2 - Generating the transfer matrix.

uresolution=0.1;

du=uresolution*1e6;
%   utargetvec=[-1:uresolution:1]*1e6;
utargetvec=[-5:uresolution:5]*1e6;


transferrow=0;
removed_rows = 0;
for phi=phivecr
    for utarget=utargetvec
        transferrow=transferrow+1;
        gamma1=acos((utarget-du/2-cos(phi/180*pi).*VPARA)./(sin(phi/180*pi).*VPERP));
        gamma2=acos((utarget+du/2-cos(phi/180*pi).*VPARA)./(sin(phi/180*pi).*VPERP));
        wv=real((gamma1-gamma2))/pi/du*dvpara*dvperp;
        transfermatrixCTS(transferrow,:)=reshape(wv,1,rows*columns);
        %if sum(sum(wv))<eps
        %    removed_rows = removed_rows + 1;
        %    transferrow=transferrow-1;
        %end
    end
end

A = transfermatrixCTS;

%% Part 3 - Generating the spectrum, S. (or b, right hand side) b = A*X + e
%[S_noisy2,e2] = Generate_noisy_spectrum(transfermatrixCTS,fvpavpe2DbiMax,noiselevel,backgroundlevel);

S0 = A*X(:);

noise_level=0.01;

background_level=1*10^0;

background_errorlevel = sqrt(background_level);
e_min_vector = ones(size(S0))*background_errorlevel;


S_noisy = S0 + noise_level*mean(sqrt(S0))*randn(size(S0)).*max([sqrt(S0) e_min_vector],[],2);

% The calculated errorbar sizes calculated from the noisy signal such that
% we can not just take the square of the errorbars and cheat that way.
e = noise_level*mean(sqrt(abs(S_noisy)))*max([sqrt(abs(S_noisy)) e_min_vector],[],2);

b = S_noisy;

%Dig into this noise. 

%% Part 4 - Solving for X.

%[xalpha,alpha] = TikhonovNonNeg(A,b,vpara,vperp,1);

%Use lsqlin instead.

%% Using MOSEK
addpath(genpath('../mosek'))

%Constructing C and d for the MOSEK routine.
%Try the different alpha parameters.
[m,n] = size(A);
alpha = 1;

X2 = X(:)/norm(X(:));

C = [A ; sqrt(alpha) * eye(n)];

%Using b = A*X doesn't help.
d = [A*X2 ; zeros(n,1)];

[x, resnorm, residual, exitflag, output] = lsqnonneg(C, d);

%figure
imagesc(reshape(x, vperpdim, vparadim)); axis xy; axis image; colorbar()

%%
figure
subplot(1,2,1)
imagesc(reshape(X, vperpdim, vparadim));
colorbar(); axis xy; axis image;

subplot(1,2,2)
imagesc(reshape(x, vperpdim, vparadim));
colorbar(); axis xy; axis image; %Good axes



