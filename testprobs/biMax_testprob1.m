%% Bi-Maxwellian.
clear, clc, close all

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

%Observation angles
phi=[10 20 40 70 85];

%Construct the u-vector.
u = construct_uvec('umax', 5*1e6, 'du', 1e5);

%Construct the normal grid and our true solution.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(vparadim, vperpdim,'vperpmin',vperpmin,'vperpmax',vperpmax,'vparamin',vparamin,'vparamax',vparamax);
xtrue = biMaxx(vpara,vperp);

%Forward model matrix A
A = transferMatrix(vpara,vperp,phi,u);

[b1, e1] = generate_noisy_b(A,xtrue);     %Basically b = A*x + noise 

vparadimfine = 100; vperpdimfine = 50;
[vparafine, vperpfine, ginfofine] = construct_vgrid(vparadimfine, vperpdimfine,'vparamax',vparamax,'vparamin',vparamin,'vperpmin',vperpmin,'vperpmax',vperpmax);
xfine = biMaxx(vparafine, vperpfine); 
Afine = transferMatrix(vparafine,vperpfine,phi,u); 
[b2, e2] = generate_noisy_b(Afine,xfine); %Same, but with a fine grid.

gu = biMaxg(u,phi);
du = u(2)-u(1);
b3 = du*gu;
[b3, e3] = add_noise(b3);

%% Plot the three right hand sides
figure
plot(u,b1(1:101))
hold on
plot(u,b2(1:101))
hold on 
plot(u,b3(1:101))
lgd = legend('Inverse crime', 'No inverse crime', 'Analytic','Interpreter','latex'); lgd.FontSize = 15;
xlabel('u','FontSize',20,'Interpreter','latex')

%%




