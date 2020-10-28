%% John vs Mosek
%Clear workspace, console etc.
clear, clc, close all

%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();
vparadim = length(gridinfo.vpara_ax);
vperpdim = length(gridinfo.vperp_ax);

%Evaluate the bi-Maxwellian on this grid with default values.
[x_true, xinfo] = biMaxx(vpara, vperp); x_true = x_true(:);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;

%ubroadening is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
ubroadening = 3; 

%Observation angles
phi=[10 20 40 70 85];

[b, binfo] = biMaxb(ustruct,phi);

u = construct_uvec(ustruct);

A = transferMatrix(vpara,vperp,phi,u);

%alpha values for 0th order Tikhonov and 1st order Tikhonov.
alphavec = logspace(5, 15, 20);

%Add noise to b.
[b_noisy, ~, e] = add_noise(b,0.01);

%This is for the 1st order Tikhonov.
L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.

alpha = alphavec(1);

%%
xJohn  = GPCG_TikhNN(A,b_noisy,alpha,L);


% Use in mosekopt
%[r,resp] = mosekopt('minimize', 'lsqnonneg', param);

%%
[m,n] = size(A);

%%
%xmosek = mosek_TikhNN(A,b_noisy,alpha,L);

param.MSK_IPAR_PRESOLVE_LINDEP_USE = 'MSK_OFF';
param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';

prob.c = [zeros(n,1); 1];
prob.a = sparse([A_IC',  ones(m-s,1);
                 A_IC', -ones(m-s,1);
                 A_I' ,  zeros(s, 1)  ]);
prob.blc = [zeros(m-s,1); -inf*ones(m-s,1); rhsI];
prob.buc = [inf*ones(m-s,1); zeros(m-s,1);  rhsI];
prob.blx = [];
prob.bux = [];

echo_level = 3;

cmd = 'minimize';
cmd = sprintf('%s echo(%i)',cmd,echo_level);

[r,res] = mosekopt(cmd,prob,param);






%% Objective function for the two

figure(1)
showDistribution(xJohn,gridinfo)

figure(2)
showDistribution(xmosek,gridinfo)

%% 
%In John's objective function, 1/2 b'*b is not in it.
f1 = @(x) 1/2*x'*(A'*A+alpha*L'*L)*x - x'*A'*b_noisy + 1/2*b_noisy'*b_noisy; %John's objective function
f2 = @(x) 1/2*norm([A;sqrt(alpha)*L]*x - [b_noisy;zeros(size(L,1),1)])^2; %Mosek's objective function

disp('Johns solution')
f1(xJohn)
f2(xJohn)

disp('moseks solution')
f1(xmosek)
f2(xmosek)

%%
clc
[M,N] = size(A);
alpha = alphavec(1);

L = reguL(vparadim,vperpdim); %L'L is eq. (16) in Jacobsen 2016 Phys Control.
%L = speye(N);
%alpha = 0;
C = [A ; sqrt(alpha)*L];
d = [b ; zeros(size(L,1),1)];

options = mskoptimset('Diagnostics','on','Display','on');

[xmosek2, resnorm, residual, exitflag, output, lambda] = lsqnonneg(C, d, zeros(size(x_true)),options);



