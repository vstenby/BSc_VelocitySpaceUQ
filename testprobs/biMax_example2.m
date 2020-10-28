%% BiMax example with three angles.
clear, clc, close all

%Load the simulation
%load('quicksim.mat');
%xsim = x; 

%Add the dependencies
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the basic grid.
vparadim = 40; vperpdim = 20;
[vpara, vperp, ginfo] = construct_vgrid(40,20);
[vparafine, vperpfine, ginfofine] = construct_vgrid(100,50);

%Construct the default u vector.
%u = construct_uvec('umin', [], 'umax', 5*1e6, 'du', 1e5);

uresolution = 0.1;
du=uresolution*1e6;
    %   utargetvec=[-1:uresolution:1]*1e6;
u=[-5:uresolution:5]*1e6;
    
    
x = biMaxx(vpara,vperp);
xfine = biMaxx(vparafine,vperpfine);

%Observation angles
%phi = [10 20 40 70 85];
%phi = [60 80 45];
phi =[10 20 40 70 85];

%Construct the A matrix.
A = transferMatrix(vpara,vperp,phi,u);
Afine = transferMatrix(vparafine,vperpfine,phi,u);

[b, e] = Generate_noisy_spectrum(Afine,xfine,0.01,1);
%[b, e] = generate_noisy_b(A,x);

%noiselevel=0.01;
%backgroundlevel=1*10^0;
%[b2, e] = Generate_noisy_spectrum(A,x(:),noiselevel,backgroundlevel)
%L = reguL(vparadim,vperpdim);


%for isigma=1:length(b)
%    A(isigma,:)=A(isigma,:)/e(isigma);
%    b(isigma)=b(isigma)/e(isigma); 
%end


%%
N = size(A,2);
%[xsim, alpha] = NNHGS(A,b,L,1000,'solver', 'lsqnonneg','welford',true,'nburnin',100);



%%
alpha = logspace(-15, -2, 30);
[~,N] = size(A);

for i=1:length(alpha)
    C = [A       ; sqrt(alpha(i))*speye(N)];
    d = [b ; zeros(N,1)];
    x0th(:,i) = lsqnonneg(C,d);
    
    C = [A       ; sqrt(alpha(i))*L];
    d = [b ; zeros(size(L,1),1)];
    x1st(:,i) = lsqnonneg(C,d);
end

%% 

figure
semilogx(alpha,relerr(x(:),x0th))
hold on
semilogx(alpha,relerr(x(:),x1st))
legend('0th order', '1st order')







