%plot a bi-Maxwellian distribution ignoring potential energy
clear;
Mp   = 1.6726e-27;
Qe   = 1.6021917e-19;
Mi=4*Mp;

%(vpara,vperp) grid
vparamax=1.4e7;
vparamin=-1.4e7;
vperpmin=1e5;
vperpmax=1.4e7;

% vparamax=4e6;
% vparamin=-4e6;
% vperpmin=1e3;
% vperpmax=4e6;

vparadim=100;
vperpdim=50;
dvpara=(vparamax-vparamin)/(vparadim-1);
dvperp=(vperpmax-vperpmin)/(vperpdim-1);

vpara=vparamin:dvpara:vparamax;
vperp=vperpmin:dvperp:vperpmax;
[VPARA,VPERP]=meshgrid(vpara,vperp);


%(E,p) grid
%boundaries of the (E,p)-space
Emin=10e3;
Emax=4e6;
%Emax=100e3;
pmin=-0.999;
pmax=0.999;
vmaxE=sqrt(2*Emax*Qe/Mi);
vminE=sqrt(2*Emin*Qe/Mi);

%resolution of the tomography grid in (E,p)
Edim=200;
pdim=200;

%coarse and fine grids for calculating synthetic data and the inversion
dE=(Emax-Emin)/(Edim-1);
dp=(pmax-pmin)/(pdim-1);

energy = Emin:dE:Emax;
pitch = pmin:dp:pmax;
[ENERGY, PITCH] = meshgrid(energy, pitch);



%plasma parameters, isotropic means Tpara=Tperp
Tpara=200e3; %eV
Tperp=200e3; %eV
ni=1e19;     %1/m^3
vparadrift=5e6;

%define thermal velocities
vthpara=sqrt(2*Tpara*Qe/Mi);
vthperp=sqrt(2*Tperp*Qe/Mi);




%3D bi-Maxwellian with drift in (vpa,vpe)
for j=1:length(vperp)
    fvpavpe3DbiMaxDrift(j,:)=ni/(pi^(3/2)*vthpara*vthperp^2)*exp(-((vpara-vparadrift)/vthpara).^2-(vperp(j)/vthperp)^2);
end

%2D bi-Maxwellian with drift in (vpa,vpe)
fvpavpe2DbiMaxDrift=fvpavpe3DbiMaxDrift.*(2*pi*VPERP);

%2D bi-Maxwellian with drift in (E,p)
fEpbiMaxDrift=ni*sqrt(ENERGY/(pi*Tpara*Tperp^2*Qe^2))...
    .*exp(-(PITCH.^2.*ENERGY*Qe+0.5*Mi*vparadrift^2 ...
    -vparadrift*PITCH.*sqrt(2*Mi*ENERGY*Qe))/(Tpara*Qe)...
    -(1-PITCH.^2).*ENERGY/Tperp)*Qe;

%slowing-down
Ecrit=44*20000;  %eV
vcrit=sqrt(2*Ecrit*Qe/Mi);
Ebirth=3.5e6; %eV
vbirth=sqrt(2*Ebirth*Qe/Mi);
Ebirthwidth=6e4;

%slowing-down in 3D and 2D (vpa,vpe)
fvpavpe3DSD=ni*3/(4*pi)/(log(1+(vbirth/vcrit)^3))./((VPARA.^2+VPERP.^2).^1.5+vcrit^3).*erfc((0.5*Mi*(VPARA.^2+VPERP.^2)/Qe-Ebirth)/Ebirthwidth)/2;
fvpavpe2DSD=fvpavpe3DSD.*(2*pi*VPERP);

%slowing down in (E,p)
fEpSD=3*ni/(4*log(1+(Ebirth/Ecrit)^1.5))*ENERGY.^0.5./(ENERGY.^1.5+Ecrit^1.5).*erfc((ENERGY-Ebirth)/Ebirthwidth)/2;
fEnergySpectrumNum=trapz(pitch,fEpSD);
fEnergySpectrum=3*ni/(2*log(1+(Ebirth/Ecrit)^1.5))*energy.^0.5./(energy.^1.5+Ecrit^1.5).*erfc((energy-Ebirth)/Ebirthwidth)/2;

if 0
%beam
Mibeam=2*Mp;
Ecritbeam=400e3;
vcritbeam=sqrt(2*Ecritbeam*Qe/Mibeam);
Ebirthbeamfull=1000e3; %eV
vbirthbeamfull=sqrt(2*Ebirthbeamfull*Qe/Mibeam);
vbirthbeamhalf=sqrt(2*Ebirthbeamfull/2*Qe/Mibeam);
vbirthbeamthird=sqrt(2*Ebirthbeamfull/3*Qe/Mibeam);
Ebirthwidthbeam=1e5;
vbirthwidthbeam=1e5;
v=vminE:(vmaxE/100):vmaxE;
[V,PITCH2]=meshgrid(v,pitch);
pbirthbeam=0.6;

VRATIOBCFULL=((V.^3/vbirthbeamfull^3).*(vbirthbeamfull^3+vcritbeam^3)./(V.^3+vcritbeam^3)).^(1/6);
VRATIOBCHALF=((V.^3/vbirthbeamhalf^3).*(vbirthbeamhalf^3+vcritbeam^3)./(V.^3+vcritbeam^3)).^(1/6);
VRATIOBCTHIRD=((V.^3/vbirthbeamthird^3).*(vbirthbeamthird^3+vcritbeam^3)./(V.^3+vcritbeam^3)).^(1/6);

%for n=0:11
SUMLEGENDREFULL=zeros(size(VRATIOBCFULL));
SUMLEGENDREHALF=zeros(size(VRATIOBCFULL));
SUMLEGENDRETHIRD=zeros(size(VRATIOBCFULL));
for n=0:20
    legendrepitch=legendreP(n,PITCH2);
    legendrepitchbirth=legendreP(n,pbirthbeam);
    SUMLEGENDREFULL=SUMLEGENDREFULL+(n+0.5)*legendrepitch.*legendrepitchbirth.*...
        VRATIOBCFULL.^(n*(n+1)).*erfc((V-vbirthbeamfull)/vbirthwidthbeam)/2;
    SUMLEGENDREHALF=SUMLEGENDREHALF+(n+0.5)*legendrepitch.*legendrepitchbirth.*...
    +VRATIOBCHALF.^(n*(n+1)).*erfc((V-vbirthbeamhalf)/vbirthwidthbeam)/2;
    SUMLEGENDRETHIRD=SUMLEGENDRETHIRD+(n+0.5)*legendrepitch.*legendrepitchbirth.*...
    +VRATIOBCTHIRD.^(n*(n+1)).*erfc((V-vbirthbeamthird)/vbirthwidthbeam)/2;
end
SUMLEGENDREFULL(find(SUMLEGENDREFULL<0))=0;
SUMLEGENDREHALF(find(SUMLEGENDREHALF<0))=0;
SUMLEGENDRETHIRD(find(SUMLEGENDRETHIRD<0))=0;
fvpbeam=1e30./(V.^3+vcritbeam^3).*(0.5*SUMLEGENDREFULL+0.3*SUMLEGENDREHALF+0.2*SUMLEGENDRETHIRD);


v_on_E_grid = sqrt(2*ENERGY*Qe/Mibeam);
fEpbeam=interp2(v,pitch,fvpbeam./(Mibeam*V),v_on_E_grid,PITCH,'linear',0);
figure(300);hold on;
%subplot(3,4,n+1)
imagesc(energy/1000,pitch,fEpbeam);colorbar; 
caxis([0 max(max(fEpbeam/10))])
set(gca,'ydir','normal')
end
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% coordinate transformations
% (E,p) -> (vpa,vpe)

%Jacobian(E,p) -> (vpa,vpe)
JacobianEpvpavpe=zeros(size(VPARA));
for jp=1:length(vperp)
    JacobianEpvpavpe(jp,:)=Mi*vperp(jp)./sqrt(vpara.^2+vperp(jp)^2)/Qe;
end

%here plus sign from Ben Geiger uses this
for jp= 1: length(vperp)
    E_on_vpavpeMESH(jp,:) =  0.5 * (vperp(jp)^2+vpara.^2)*Mi/Qe;
    p_on_vpavpeMESH(jp,:) = vpara./sqrt(vperp(jp)^2+vpara.^2);
end;

%interpolate (E,p) -> (vpa,vpe)
fvpavpe2DtransformedbiMax=interp2(energy,pitch,fEpbiMaxDrift,E_on_vpavpeMESH,p_on_vpavpeMESH,'*spline');
%fvpavpe2DtransformedbiMax(fvpavpe2DtransformedbiMax < 1e-6) = 0;
fvpavpe2DtransformedbiMax=fvpavpe2DtransformedbiMax.*JacobianEpvpavpe;
fvpavpe3DtransformedbiMax=fvpavpe2DtransformedbiMax./(2*pi*VPERP);

fvpavpe2DtransformedSD=interp2(energy,pitch,fEpSD,E_on_vpavpeMESH,p_on_vpavpeMESH,'*spline');
%fvpavpe2DtransformedSD(fvpavpe2DtransformedSD < 1e-6) = 0;
fvpavpe2DtransformedSD=fvpavpe2DtransformedSD.*JacobianEpvpavpe;
fvpavpe3DtransformedSD=fvpavpe2DtransformedSD./(2*pi*VPERP);

% (vpa,vpe) -> (E,p)
JacobianvpavpeEp = 1/Mi*sqrt(VPARA.^2+VPERP.^2)./VPERP*Qe;

v_par_on_E_grid = PITCH.*sqrt(2*ENERGY*Qe/Mi);
v_perp_on_E_grid = sqrt(1-PITCH.^2).*sqrt(2*ENERGY*Qe/Mi);

fEptransformedbiMax= interp2(vpara,vperp,fvpavpe2DbiMaxDrift.*JacobianvpavpeEp,v_par_on_E_grid,v_perp_on_E_grid,'linear',0);
%fEptransformedbiMax(fEptransformedbiMax < 1e-6) = 0;

fEptransformedSD= interp2(vpara,vperp,fvpavpe2DSD.*JacobianvpavpeEp,v_par_on_E_grid,v_perp_on_E_grid,'linear',0);
%fEptransformedSD(fEptransformedSD < 1e-6) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate densities and differences
fEpdifferencebiMax=fEpbiMaxDrift-fEptransformedbiMax;
fvpavpe2DdifferencebiMax=fvpavpe2DbiMaxDrift-fvpavpe2DtransformedbiMax;
fvpavpe3DdifferencebiMax=fvpavpe3DbiMaxDrift-fvpavpe3DtransformedbiMax;

fEpdifferenceSD=fEpSD-fEptransformedSD;
fvpavpe2DdifferenceSD=fvpavpe2DSD-fvpavpe2DtransformedSD;
fvpavpe3DdifferenceSD=fvpavpe3DSD-fvpavpe3DtransformedSD;

%integrate to get fast-ion density
nbiMaxDriftvpavpe=trapz(vperp',trapz(vpara',fvpavpe2DbiMaxDrift'));
nbiMaxDriftEp = trapz(pitch,trapz(energy,fEpbiMaxDrift'));
nbiMaxDriftTransformedvpavpe = trapz(vperp,trapz(vpara,fvpavpe2DtransformedbiMax'));
nbiMaxDriftTransformedEp = trapz(pitch,trapz(energy,fEptransformedbiMax'));
nbiMaxDriftdifferencevpavpe=trapz(vperp,trapz(vpara,fvpavpe2DdifferencebiMax'));
nbiMaxDriftdifferenceEp=trapz(pitch,trapz(energy,fEpdifferencebiMax'));

nSDvpavpe=trapz(vperp',trapz(vpara',fvpavpe2DSD'));
nSDEp = trapz(pitch,trapz(energy,fEpSD'));
nSDTransformedvpavpe = trapz(vperp,trapz(vpara,fvpavpe2DtransformedSD'));
nSDTransformedEp = trapz(pitch,trapz(energy,fEptransformedSD'));
nSDdifferencevpavpe=trapz(vperp,trapz(vpara,fvpavpe2DdifferenceSD'));
nSDdifferenceEp=trapz(pitch,trapz(energy,fEpdifferenceSD'));
nSDEnergySpectrum=trapz(energy,fEnergySpectrum);
nSDEnergySpectrumNum=trapz(energy,fEnergySpectrumNum);

%observation angles
phivec=[40];

%number of points per spectrum
udim=200;

%spectral resolution of the measurements divided by bin width u of the
%spectra
ubroadening=3;

%spectral data points
umin=sqrt(2*Emin*Qe/Mi);
umax=sqrt(2*Emax*Qe/Mi);
du=(umax-umin)/(udim/2-1);
dures=ubroadening*du;
%uvec=[-umax:du:-umin umin:du:umax];
uvec=[-umax:du:umax];

Wrow=0;
disp('computing transfer matrix...');
for iphi=1:length(phivec);
    phi=phivec(iphi);
    for iu=1:length(uvec)
        u=uvec(iu);
        Wrow=Wrow+1;
        gamma1=acos((u-dures./2-cos(phi/180*pi).*VPARA)./(sin(phi/180*pi).*VPERP));
        gamma2=acos((u+dures./2-cos(phi/180*pi).*VPARA)./(sin(phi/180*pi).*VPERP));
        wv=real((gamma1-gamma2))/pi/dures*dvpara*dvperp;
        Wvpavpe(Wrow,:)=reshape(wv,1,vparadim*vperpdim);
        if sum(sum(wv))<eps
            %Wrow=Wrow-1;
        end
    end
    
end
phivec=[40];
Wrow=0;
for iphi=1:length(phivec);
    phi=phivec(iphi);
    for iu=1:length(uvec)
        u=uvec(iu);
        Wrow=Wrow+1;
        gamma1Ep = acos(1./sin(phi.*pi./180).*1./sqrt(1-PITCH.^2).*((u-dures./2)./sqrt(2.*ENERGY./Mi*Qe)-PITCH.*cos(phi.*pi./180)));
        gamma2Ep = acos(1./sin(phi.*pi./180).*1./sqrt(1-PITCH.^2).*((u+dures./2)./sqrt(2.*ENERGY./Mi*Qe)-PITCH.*cos(phi.*pi./180)));
        wEp=real((gamma1Ep-gamma2Ep))/pi/dures*dE*dp;
        WEp(Wrow,:)=reshape(wEp,1,Edim*pdim);
        if sum(sum(wEp))<eps
            %Wrow=Wrow-1;
        end
    end
    
end



Teff=Tperp*(sin(phi/180*pi))^2+Tpara*(cos(phi/180*pi))^2;

gubiMaxanalytic=ni*(Mi/(2*pi*Teff*Qe))^0.5*exp(-(Mi*(uvec-vparadrift*cos(phi/180*pi)).^2)/(2*Teff*Qe));
gubiMaxnumvpavpe=Wvpavpe*reshape(fvpavpe2DbiMaxDrift,vparadim*vperpdim,1);

%% Per Christian plot

plot(gubiMaxanalytic)
hold on
plot(gubiMaxnumvpavpe)


%%
gubiMaxnumvpavpetransformed=Wvpavpe*reshape(fvpavpe2DtransformedbiMax,vparadim*vperpdim,1);
gubiMaxnumEp=WEp*reshape(fEpbiMaxDrift,Edim*pdim,1);
gubiMaxnumEptransformed=WEp*reshape(fEptransformedbiMax,Edim*pdim,1);


%gunaivSD=0.5*ni/(log((vbirth+vcrit)^2/(vbirth^2+vcrit^2-vbirth*vcrit))-sqrt(12)*atan((1-2*vbirth/vcrit)/sqrt(3))+pi/sqrt(3))...
%            *3*vcrit^2./(abs(uvec).^3+vcrit^3).*erfc((0.5*Mi*uvec.^2/Qe-Ebirth)/Ebirthwidth); %0.5 in the beginning due to integration from 0 to infty
 
guSD=ni/((log(1+(vbirth/vcrit)^3))*4*vcrit)*(log(abs((vbirth^2-vcrit*vbirth+vcrit^2)./(uvec.^2-vcrit*abs(uvec)+vcrit^2)))+...
    log(((abs(uvec)+vcrit)/(vbirth+vcrit)).^2)+2*sqrt(3)*(atan((2*vbirth-vcrit)/(sqrt(3)*vcrit))-atan((2*abs(uvec)-vcrit)/(sqrt(3)*vcrit)))).*erfc((0.5*Mi*uvec.^2/Qe-Ebirth)/Ebirthwidth)/2; 
guSDnumvpavpe=Wvpavpe*reshape(fvpavpe2DSD,vparadim*vperpdim,1);
guSDnumvpavpetransformed=Wvpavpe*reshape(fvpavpe2DtransformedSD,vparadim*vperpdim,1);
guSDnumEp=WEp*reshape(fEpSD,Edim*pdim,1);
guSDnumEptransformed=WEp*reshape(fEptransformedSD,Edim*pdim,1);


vmag=[du:du:umax];
fvSDAnalytic=3*ni/(4*pi*log(1+(vbirth/vcrit)^3)) ...
       ./(abs(vmag).^3+vcrit^3).*erfc((0.5*Mi*vmag.^2/Qe-Ebirth)/Ebirthwidth)/2;
guSDPlus=ni/((log(1+(vbirth/vcrit)^3))*4*vcrit)*(log(abs((vbirth^2-vcrit*vbirth+vcrit^2)./(vmag.^2-vcrit*abs(vmag)+vcrit^2)))+...
    log(((abs(vmag)+vcrit)/(vbirth+vcrit)).^2)+2*sqrt(3)*(atan((2*vbirth-vcrit)/(sqrt(3)*vcrit))-atan((2*abs(vmag)-vcrit)/(sqrt(3)*vcrit)))).*erfc((0.5*Mi*vmag.^2/Qe-Ebirth)/Ebirthwidth)/2; 

A=[];

for i=1:length(vmag)
    A=[A;vmag*2*pi*du];
end
for i=2:length(vmag)
    A(i,1:i-1)=0;
end
guSDPlusnumfv1D=A*fvSDAnalytic';


nSDgu=trapz(uvec',guSD');
nSDgunumvpavpe=trapz(uvec',guSDnumvpavpe');
nSDgunumvpavpetransformed=trapz(uvec',guSDnumvpavpetransformed');
nSDgunumEp=trapz(uvec',guSDnumEp');
nSDgunumEptransformed=trapz(uvec',guSDnumEptransformed');
nSDguPlus=2*trapz(vmag',guSDPlus');
nSDguPlusnumfv1D=2*trapz(vmag',guSDPlusnumfv1D');




Evmag=1/2*Mi*vmag.^2/Qe;


figure(10);clf;hold on;
%plot(uvec,guSD);
subplot(1,3,1)
hold on;box on;
plot(vmag,guSDPlus);
%plot(vmag,guSDPlusnumfv1D);
%plot(vmag,guSDPlusNoise);
% plot(vmag,guSDPlusnumfv1DNoise);
%plot(-vmag,guSDPlusnumfv1D);

subplot(1,3,2)
hold on;box on;
plot(vmag,fvSDAnalytic);

subplot(1,3,3)
hold on;box on;
plot(energy,fEnergySpectrum);



if 1
figure(20);clf;hold on;
%plot(uvec,gunaivSD,'LineWidth',2)
plot(uvec,guSD,'LineWidth',3)
plot(uvec,guSDnumvpavpe,'LineWidth',2)
plot(uvec,guSDnumEp,'--','LineWidth',2)
plot(vmag,guSDPlusnumfv1D,':')
plot(-vmag,guSDPlusnumfv1D,':')
% plot(uvec,gunumSDvpavpetransformed,'s')
% plot(uvec,gunumSDEptransformed,'+')
%plot(uvec(end/2+1:end),G1DEnergySpectrum,'g')


figure(30);clf;hold on;
plot(uvec,gubiMaxanalytic,'LineWidth',2)
plot(uvec,gubiMaxnumvpavpe,'o')
plot(uvec,gubiMaxnumEp,'x')
plot(uvec,gubiMaxnumvpavpetransformed,'s')
plot(uvec,gubiMaxnumEptransformed,'+')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots

fs=12;
axvpavpe = [-14 14 -0.05 14];
axEp = [-0.01 4 -1.01 1.01];
isocontours=[1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10];

load MyColormap.mat



if 1
fig=figure(100);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
colormap(mycolormap)
subplot(4,3,1)
%contour(vpara/1e6,vperp/1e6,fvpavpe3DbiMaxDrift,isocontours);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe3DbiMaxDrift,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^3/m^6]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe3DbiMaxDrift))])

subplot(4,3,4)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe2DbiMaxDrift))])

subplot(4,3,7)
%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,fEpbiMaxDrift,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fEpbiMaxDrift))])
% set(handlecolorbar,'fontsize',fs,'fontweight','bold');
% set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);


subplot(4,3,2)
%contour(vpara/1e6,vperp/1e6,fvpavpe3DbiMaxDrift,isocontours);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe3DtransformedbiMax,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^3/m^6]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe3DbiMaxDrift))])


subplot(4,3,5)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe2DtransformedbiMax,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe2DbiMaxDrift))])

subplot(4,3,8)
%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,fEptransformedbiMax,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fEpbiMaxDrift))])


subplot(4,3,6)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe2DdifferencebiMax,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
%caxis([-max(max(fvpavpe2DbiMaxDrift))/100 max(max(fvpavpe2DbiMaxDrift))/100])


subplot(4,3,9)
%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,fEpdifferencebiMax,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
%caxis([-max(max(fEpbiMaxDrift))/100 max(max(fEpbiMaxDrift))/100])
% set(handlecolorbar,'fontsize',fs,'fontweight','bold');
% set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
end


if 1
fig=figure(200);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
colormap(mycolormap)

subplot(4,3,1)
%contour(vpara/1e6,vperp/1e6,fvpavpe3DbiMaxDrift,isocontours);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe3DSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^3/m^6]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe3DSD))])

subplot(4,3,4)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe2DSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe2DSD))])

subplot(4,3,7)
%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,fEpSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fEpSD))])
% set(handlecolorbar,'fontsize',fs,'fontweight','bold');
% set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);


subplot(4,3,2)
%contour(vpara/1e6,vperp/1e6,fvpavpe3DbiMaxDrift,isocontours);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe3DtransformedSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^3/m^6]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
caxis([0 max(max(fvpavpe3DSD))])
handlecolorbar=colorbar;%('horiz');

subplot(4,3,5)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,fvpavpe2DtransformedSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fvpavpe2DSD))])

subplot(4,3,8)
%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,fEptransformedSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fEpSD))])

subplot(4,3,3)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,abs(fvpavpe3DdifferenceSD),50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
%caxis([-max(max(fvpavpe3DSD))/100 max(max(fvpavpe3DSD))/100])
caxis([0 max(max(fvpavpe3DSD))/1e3])

subplot(4,3,6)
%contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMaxDrift,isocontours*1e6);
[dummy,H]=contourf(vpara/1e6,vperp/1e6,abs(fvpavpe2DdifferenceSD),50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(axvpavpe)
text(1, 2.6,'f [s^2/m^5]','fontsize',fs,'fontweight','bold')
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
%caxis([-max(max(fvpavpe2DSD))/100 max(max(fvpavpe2DSD))/100])
caxis([0 max(max(fvpavpe2DSD))/1e3])

subplot(4,3,9)
%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,abs(fEpdifferenceSD),50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
%caxis([-max(max(fEpSD))/100 max(max(fEpSD))/100])
caxis([0 max(max(fEpSD))/1e3])
% set(handlecolorbar,'fontsize',fs,'fontweight','bold');
% set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
end

if 0
%Alfven resonances
B=5.3;
mu0=4*pi*1e-7;
n=1e20;
Mi=2.5*1.67e-27;

vA=B/sqrt(mu0*n*Mi);

pAlfven=-1:0.01:1
EAlfven1=0.5*2*1.67e-27*vA^2./(pAlfven.^2)/1.6e-19/1e6;
EAlfven2=0.5*2*1.67e-27*(vA/3)^2./(pAlfven.^2)/1.6e-19/1e6;
    
fig=figure(210);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
colormap(mycolormap)

%contour(E_coarse/1000,p_coarse,fEpbiMaxDrift,isocontours*1e14);
[dummy,H]=contourf(energy/1e6,pitch,fEpSD,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
plot(EAlfven1,pAlfven,'c')
plot(EAlfven2,pAlfven,'m')
axis('square')
axis(axEp)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [MeV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
handlecolorbar=colorbar;%('horiz');
caxis([0 max(max(fEpSD))])
% set(handlecolorbar,'fontsize',fs,'fontweight','bold');
% set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
end