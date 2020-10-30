%plot a bi-Maxwellian distribution ignoring potential energy
clear;

Mp   = 1.6726e-27;
Qe   = 1.6021917e-19;
MD=2*Mp;

%plasma parameters, isotropic means Tpara=Tperp
Tpara=20; %keV
Tperp=20; %keV
ne=1e19;     %1/m^3
vdrift=5e5; %m/s

%define thermal velocities
vthpara=sqrt(2*Tpara*Qe*1000/MD);
vthperp=sqrt(2*Tperp*Qe*1000/MD);

%set the grids. dvfine for computation of spectra, dv for inversion
%(vpara,vperp) grid
vparamax=4e6;
vparamin=-4e6;
vperpmin=1e4;
vperpmax=4e6;


vparadim=40;
vperpdim=20;
dvpara=(vparamax-vparamin)/(vparadim-1);
dvperp=(vperpmax-vperpmin)/(vperpdim-1);

vpara=vparamin:dvpara:vparamax;
vperp=vperpmin:dvperp:vperpmax;
[VPARA,VPERP]=meshgrid(vpara,vperp);
[rows,columns]=size(VPERP);

%(vparafine,vperpfine) grid
vparadimfine=100;
vperpdimfine=50;
dvparafine=(vparamax-vparamin)/(vparadimfine-1);
dvperpfine=(vperpmax-vperpmin)/(vperpdimfine-1);

vparafine=vparamin:dvparafine:vparamax;
vperpfine=vperpmin:dvperpfine:vperpmax;
[VPARAFINE,VPERPFINE]=meshgrid(vparafine,vperpfine);
[rowsfine,columnsfine]=size(VPERPFINE);



%bi-Maxwellian distribution function as slice in 3D
fvpavpe3DbiMaxfine=ne/(pi^(3/2)*vthpara*vthperp^2)*exp(-((VPARAFINE-vdrift)/vthpara).^2-(VPERPFINE/vthperp).^2);
fvpavpe3DbiMax=ne/(pi^(3/2)*vthpara*vthperp^2)*exp(-((VPARA-vdrift)/vthpara).^2-(VPERP/vthperp).^2);

%fvpavpe in 2D
fvpavpe2DbiMax=fvpavpe3DbiMax.*(2*pi*VPERP);
fvpavpe2DbiMaxfine=fvpavpe3DbiMaxfine.*(2*pi*VPERPFINE);

%integrate to get fast-ion density
%nvpavpebiMax=trapz(vperpfine',trapz(vparafine',fvpavpe2DbiMax'));

%observation angles
phivecr=[10 20 40 70 85];

%This is equivalent to figure 4 in Bi-Maxwellian document.
%imagesc(vpara/1e6,vperp/1e6,reshape(fvpavpe2DbiMax,length(vperp),length(vpara)))

for uresolution=0.1
    disp(uresolution)
    du=uresolution*1e6;
    %   utargetvec=[-1:uresolution:1]*1e6;
    utargetvec=[-5:uresolution:5]*1e6;
    disp(du)
    %   for phivecr=2:2:90
    %for phi1=4:4:90
    %    for phi2=4:4:phi1
    %phivecr=[phi1 phi2];
    transferrow=0;
    disp('computing transfer matrix...');
    if 1
        for phi=phivecr;
            for utarget=utargetvec
                transferrow=transferrow+1;
                gamma1fine=acos((utarget-du/2-cos(phi/180*pi).*VPARAFINE)./(sin(phi/180*pi).*VPERPFINE));
                gamma2fine=acos((utarget+du/2-cos(phi/180*pi).*VPARAFINE)./(sin(phi/180*pi).*VPERPFINE));
                wvfine=real((gamma1fine-gamma2fine))/pi/du*dvparafine*dvperpfine;
                transfermatrixCTSfine(transferrow,:)=reshape(wvfine,1,rowsfine*columnsfine);
                gamma1=acos((utarget-du/2-cos(phi/180*pi).*VPARA)./(sin(phi/180*pi).*VPERP));
                gamma2=acos((utarget+du/2-cos(phi/180*pi).*VPARA)./(sin(phi/180*pi).*VPERP));
                wv=real((gamma1-gamma2))/pi/du*dvpara*dvperp;
                transfermatrixCTS(transferrow,:)=reshape(wv,1,rows*columns);
                %if sum(sum(wv))<eps
                %    transferrow=transferrow-1; 
                %    %Viktor: I commented this out - no idea why this is here...
                %end
            end
        end
    end
end

A = transfermatrixCTSfine;


%plot examples lines from the fine transfer matrix
%for i=10:10:375
%figure(100+i);clf;hold on;box on;set(gca, 'Layer','top');
%imagesc(vparafine/1e6,vperpfine/1e6,reshape(transfermatrixCTSfine(i,:),length(vperpfine),length(vparafine)))
%set(gca,'ydir','normal')
%colorbar
%caxis([0 max(max(transfermatrixCTSfine(i,:)))])
%axis equal
%set(gca,'Fontsize',16)
%xlabel('v_{||} [10^6 m/s]')
%ylabel('v_\perp [10^6 m/s]')
%axis([vparafine(1) vparafine(end) 0 vperpfine(end)]/1e6)
%colormap(flipud(hot))
%end

%plot examples lines from the coarse transfer matrix
% for i=10:10:375
% figure(100+i);clf;hold on;box on;set(gca, 'Layer','top');
% imagesc(vpara/1e6,vperp/1e6,reshape(transfermatrixCTS(i,:),length(vperp),length(vpara)))
% set(gca,'ydir','normal')
% colorbar
% caxis([0 max(max(transfermatrixCTS(i,:)))])
% axis equal
% set(gca,'Fontsize',16)
% xlabel('v_{||} [10^6 m/s]')
% ylabel('v_\perp [10^6 m/s]')
% axis([vpara(1) vpara(end) 0 vperp(end)]/1e6)
% colormap(flipud(hot))
% end


noiselevel=0.01;
backgroundlevel=1*10^0;



[S_noisy,e] = Generate_noisy_spectrum(transfermatrixCTSfine,fvpavpe2DbiMaxfine,noiselevel,backgroundlevel);


% fvpavpe2DbiMaxcolumn = reshape(fvpavpe2DbiMax,size(transfermatrixCTS,2),1);
% %S_noisy=transfermatrixCTS*fvpavpe2DbiMaxcolumn;
% figure;plot(S_noisy)
%e=S_noisy*1e-5+max(S_noisy)*0.001;


for isigma=1:length(S_noisy)
    transfermatrixCTS(isigma,:)=transfermatrixCTS(isigma,:)/e(isigma);
    disp('Error')
    disp(e(1))
    S_noisy_norm(isigma)=S_noisy(isigma)/e(isigma); 
end
S_noisy_norm=S_noisy_norm';


clc
disp(vparadim)
disp(vperpdim)
disp(size(fvpavpe2DbiMax))

%L = reguL2(n1,n2) d

%[xalpha,alpha] = TikhonovNonNeg(transfermatrixCTS,S_noisy_norm,vpara,vperp,1);


%%




%%
figure(10); clf; hold on; box on;set(gca, 'Layer','top');
for i=1:length(alpha)
    subplot(2,3,i)
    imagesc(vpara/1e6,vperp/1e6,reshape(xalpha(:,i),length(vperp),length(vpara)))
    set(gca,'ydir','normal')
    colorbar
%    caxis([0 max(max(fvpavpe2DbiMax))])
    axis equal
    set(gca,'Fontsize',16)
    xlabel('v_{||} [10^6 m/s]')
    ylabel('v_\perp [10^6 m/s]')
    axis([vpara(1) vpara(end) 0 vperp(end)]/1e6)
end
%set(gcf,'colormap',colormap_rel_diff)
colormap(flipud(hot))


figure(20);clf;hold on;box on;set(gca, 'Layer','top');
imagesc(vpara/1e6,vperp/1e6,reshape(fvpavpe2DbiMax,length(vperp),length(vpara)))
set(gca,'ydir','normal')
colorbar
caxis([0 max(max(fvpavpe2DbiMax))])
axis equal
set(gca,'Fontsize',16)
xlabel('v_{||} [10^6 m/s]')
ylabel('v_\perp [10^6 m/s]')
axis([vpara(1) vpara(end) 0 vperp(end)]/1e6)
colormap(flipud(hot))
title('SIMULATION')


