%% Drifting bi-Maxwellian UQ
%Commented out so it doesn't reset on HPC.
%clear, clc, close all

disp('Starting script')

%Add dependencies.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

%method = 'nonnegative'; %unconstrained or nonnegative
%Lidentity = true; %if L is the identity.

%Construct the default grid
[vpara, vperp, gridinfo] = construct_vgrid();

%Evaluate the bi-Maxwellian on this grid with default values.
[X, biMaxXinfo] = biMaxX(vpara, vperp);

%Construct the analytic projection.
%Boundaries of the (E,p)-space
ustruct.Emin = 10e3;
ustruct.Emax = 4e6;
%Number of points per spectrum
ustruct.udim = 200;
%Observation angles
phi=[10 20 40 70 85];

[S, biMaxSinfo] = biMaxS(ustruct,phi);
[S_noisy, s] = add_noise(S,0.01);

S = S_noisy;
%Generate A from vpara, vperp, u and phi from biMaxX and biMaxS.
%First argument is ubroadening, which is the spectral resolution 
%of the measurements divided by bin width u of the spectra.
A = biMaxA(3,biMaxXinfo,biMaxSinfo);

%% Normalize the problem.
%[Shat, Ahat] = measurement_normalization(S,A,s);
%[S2hat, W2hat, factor] = numeric_normalization(Shat, Ahat);

%Overwrite S and A with the new right hand side and left hand side.
%A = W2hat;
%S = S2hat;
%alpha = 3.4e-16;

alpha = 2.8e8;
xalpha = mosek_TikhNN(A, S, alpha);

a_params.alpha      = alpha;
params.a_params     = a_params;

showDistribution(xalpha, gridinfo)

[M, N] = size(A);

%Initialize lambda0 and delta0.

lambda_est = 1/norm(S(:)-A*xalpha(:))^2;
delta_est  = alpha*lambda_est;

% Initialization before sampling.  
nsamps     = 1000;
lamsamp    = zeros(nsamps,1); lamsamp(1) = lambda_est;
delsamp    = zeros(nsamps,1); delsamp(1) = delta_est;
xsamp      = zeros(N,nsamps);
xtemp      = xalpha; 
xsamp(:,1) = xtemp(:);

%Hyperpriors
a0         = 1; 
t0         = 0.0001; 
a1         = 1; 
t1         = 0.0001;


%L as the identity.
L = eye(N);


params.max_cg_iter  = 2000;
params.cg_step_tol  = 1e-6;
params.grad_tol     = 1e-6;
params.cg_io_flag   = 0;
params.cg_figure_no = [];
params.precond      = [];

Bmult               = 'Bmult_TomographyGMRF';
a_params.A          = A;
a_params.L          = L;


%% The UQ part.

lambda = lamsamp(1);
delta = delsamp(1);
disp('Starting sampler')
AtA = A'*A;
Atb = A'*S;
for i=1:nsamps-1
    disp(i)
    %------------------------------------------------------------------
    Axtemp       = A*xtemp;
    lamsamp(i+1) = gamrnd(a0+M/2,1./(t0+norm(Axtemp(:)-S(:))^2/2));
    %------------------------------------------------------------------
    % 1b. Using conjugacy, sample regularization precisions delta,
    % conjugate prior: delta~Gamma(a1,1/t1);
    Lxtemp       = xtemp;
    delsamp(i+1) = gamrnd(a1+N/2,1./(t1+xtemp(:)'*Lxtemp(:)/2));
    %------------------------------------------------------------------
    % 2. Use CG to sample x.
    %alpha = delsamp(i+1)/lamsamp(i+1);
    a_params.alpha = delsamp(i+1)/lamsamp(i+1);
    params.a_params = a_params;
    
    %For nonneg, we need to construct the rhs. of (5.12)
    eta = sqrt(lamsamp(i+1))*A'*randn(M,1) + sqrt(delsamp(i+1)).*randn(N,1);
    c   = Atb + eta/lamsamp(i+1); %RHS of (5.12)

    coefstruct.Afun = 'Bmult_TomographyGMRF';
    coefstruct.Aparams = a_params;

    %Now, we pass this to GPCG. We cannot necessarily do the trick from
    %(5.45), since x has to be greater than 0, not P. Therefore, we need to
    %use (5.12) with zero vector as initial guess.

    %For the nonneg, we don't use the parameters for the CG optimizer.
    %These parameters are taken from TwoDBlurSampleConjNonneg.
    [xtemp, ~, iter] = GPCG(coefstruct, c, zeros(N,1));
    %[xtemp, ~, iter] = GPCG(coefstruct, c, zeros(N,1), 50, 5, 20, 1e-6);
    
    
    xsamp(:,i+1) = xtemp;
    %showDistribution(xtemp, gridinfo)
    %drawnow()
end   


%% Show the chains



