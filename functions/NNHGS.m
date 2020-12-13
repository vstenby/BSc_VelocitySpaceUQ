function [x, alpha, delta, lambda, info] = NNHGS(A,b,L,n,varargin)
% Nonnegative Hierachical Gibbs Sampler
%
% Usage: 
%    ``[x, alpha] = NNHGS(A,b,L,n)`` 
%
%    ``[x, alpha, delta, lambda, info] = NNHGS(A,b,L,n,varargin)`` 
%
% Inputs:
%    * **A**:                   The system matrix A.
%
%    * **b**:                   The right hand side.
%
%    * **L**:                   The regularization matrix. ``[]`` is interpreted as ``speye(N)``
%
%    * **n**:                   Number of samples to be generated.
%
% Optional inputs:
%    * **disp_waitbar**:        Whether or not a waitbar should be displayed.
%   
%    * **welfordsalgorithm**:   If ``true``, then Welford's Online Algorithm is used to return mean and standard deviation of the samples.
%   
%    * **nburnin**:             Samples to be generated before Welford's Algorithm starts. Only used if ``welford=true``. 
%   
%    * **keep**:                Which samples we should keep. Only used if ``welford=true``. Default value is 1, meaning that every sample is kept.
%
%    * **solver**:              Which solver is used. ``lsqnonneg`` is default, but ``GPCG`` or ``\`` can be set.
%
%    * **scaling**:             Whether or not A should be rescaled such that the largest element in A = 1. Default is ``true``.
%
% Output:
%    * **x**:               	Samples drawn from the posterior
%
%    * **alpha**:               alpha-chain sampled from the posterior.
%
%    * **delta**:               delta-chain sampled from the posterior.
%
%    * **lambda**:              lambda-chain sampled from the posterior.
%
%    * **info**:                MATLAB struct containing all sorts of information.

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
%Set the default parameters for the NNHGS
disp_waitbar = true;
welfordsalgorithm = false;
keep = 1;
solver  = 'lsqlin';
scaling = true; %This should then be fixed.

%Unpack the varargin and evaluate.
validvars = {'disp_waitbar','welfordsalgorithm','keep','solver','scaling','nburnin'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end

if ~islogical(welfordsalgorithm), error('Welford should be logical'), end
if ~welfordsalgorithm && keep ~= 1
   warning('Trimming has to be done manually if Welford is false.') 
end

%if ~islogical(scaling), error('Scaling should be logical'), end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

if ~check_mosek() && strcmpi(solver,'lsqnonneg')
   warning('mosek is not installed, switching default solver to lsqlin.') 
   solver = 'lsqlin';
end

if isstruct(L) && ~strcmpi(scaling,'norm')
    scaling = 'norm';
    %If L is a struct containing L1p, L1E and L0b, then 
    %we should Birgitte's scaling method.
    %warning('L is a struct, changing scaling type to norm.')
end

%Scaling is set to 1/max(A) if need be, otherwise 1.
if strcmpi(scaling,'1/maxA')
    scaling_factor = 1/max(A(:));
    %Scale A.
    A = A*scaling_factor;
elseif strcmpi(scaling,'norm')
    [A, b, L, scaling_factor] = norm_normalization(A,b,L);
else
    scaling_factor = 1;
end

if welfordsalgorithm
    if ~exist('nburnin','var')
        error("Burn-in should be specified with Welford's Algorithm")
    end
    %Calculate the number of 
    nidx = [1:keep:n*keep];
    nidx = nidx + nburnin;
    nsamps = nidx(end);
else
    if exist('nburnin','var'), warning('nburnin specified but not used.'); end
    nidx = 1:n;
    nsamps = n;
    nburnin = 0;
end

%Allocate the delta, lambda and alpha
del_sim  = zeros(nsamps,1);
lam_sim  = zeros(nsamps,1);
alph_sim  = zeros(nsamps,1);

[M,N] = size(A);

%Lower bound for lsqlin.
lb = zeros(N,1);

if strcmpi(solver,'lsqlin')
    opts = optimset('display','off'); 
end


if all(size(L) == 0)
   L = speye(N); 
end

if welfordsalgorithm
    x = zeros(N,2);
else
    x = zeros(N,nsamps);
end

if disp_waitbar
    f = waitbar(1/nsamps,'Finding initial xalpha');
end

%Calculate L'*L beforehand.
LtL = L'*L;

%Initial sample.
alpha0 = 0; %This should be a parameter in the options as well.

xtemp = TikhNN(A,b,alpha0,L,'solver',solver);

alph_temp = alpha0;
lam_temp  = 1/norm(b(:)-A*xtemp(:))^2;
del_temp  = alpha0*lam_temp;

%Set the initial values for alpha, x, delta and lambda.
alph_sim(1)    = alph_temp;
lam_sim(1)     = lam_temp;
del_sim(1)     = del_temp;

%Hyperpriors
a0 = 1; 
t0 = 0.0001; 
a1 = 1;
t1 = 0.0001;

%Index for Welford's Online Algorithm
if welfordsalgorithm, count = 0; end

for i=2:nsamps
   if disp_waitbar, waitbar(i/nsamps, f, 'Sampling...'), end
   
   %NNHGS loop:
   Axtemp = A*xtemp;
   LtLxtemp = LtL*xtemp;
   
   %Note here that 1./ is because of the way MATLAB does gamrnd. 
   %Equation (58) and (59) in Dental CT
   lam_temp  = gamrnd(a0 + M/2, 1./(t0+norm(Axtemp(:)-b(:))^2/2)); 
   del_temp  = gamrnd(a1 + N/2, 1./(t1+xtemp(:)'*LtLxtemp(:)/2)); 
   
   alph_temp = del_temp/lam_temp;
   
   %Multivariate normal distribution vectors bhat and chat are generated.
   %https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Computational_methods
   bhat = b + sqrt(lam_temp^(-1))*speye(M)*randn(M,1);
   chat =     sqrt(del_temp^(-1))*speye(size(L,1))*randn(size(L,1),1);
      
   switch solver
       case '\'
            B = [sqrt(lam_temp)*A ; sqrt(del_temp)*L];
            d = [sqrt(lam_temp)*bhat ; sqrt(del_temp)*chat];
            xtemp = B\d;
        case 'lsqnonneg'
            B = [sqrt(lam_temp)*A ; sqrt(del_temp)*L];
            d = [sqrt(lam_temp)*bhat ; sqrt(del_temp)*chat];
            xtemp = lsqnonneg(B,d);
       case 'lsqlin'
            B = [sqrt(lam_temp)*A ; sqrt(del_temp)*L];
            d = [sqrt(lam_temp)*bhat ; sqrt(del_temp)*chat];
            xtemp = lsqlin(B,d,[],[],[],[],lb,[],zeros(N,1),opts);
        %case 'GPCG'
            %Right hand side of (2.5) in BaHa20 divided with lambda.
        %    rhs = A'*bhat + alph_temp*L'*chat;
        %    B = @(x) (A'*(A*x)) + alph_temp*LtL*x; 
        %    xtemp = GPCG(B, rhs, x0, 50, 5, 20, 1e-6);
       otherwise
         error('Wrong solver specified.')
   end

   if welfordsalgorithm
       if (mod(i-nburnin,keep) || keep == 1) == 1 && nburnin < i
           %True if i = (nburnin+1)+keep*p, where p = 0, 1, ...
           if count == 0
               Mk = xtemp;
               Sk = zeros(size(xtemp));
               count = 1;
           else
               [Mk, Sk, count] = welford(Mk, Sk, count, 'update', xtemp);
           end
       end
   else
       %If no Welford, then we should save all samples
       x(:,i) = xtemp;
   end
   %alpha, delta and lambda chains should be saved either way.
   alph_sim(i) = alph_temp;
   del_sim(i) = lam_temp;
   lam_sim(i) = del_temp;
end

%Saving to the info struct.
info.welford = welfordsalgorithm;
info.n = n;
info.nsamps = nsamps;
info.nburnin = nburnin;
info.keep = keep;
info.nidx = nidx';
info.alpha = alph_sim;
info.delta = del_sim;
info.lambda = lam_sim;

if welfordsalgorithm   
    %Remove the burnin and trim the three chains.
    alpha = alph_sim(info.nidx);
    delta = del_sim(info.nidx);
    lambda = lam_sim(info.nidx);
    [xmean, xstd] = welford(Mk, Sk, count, 'finalize');
    x = [xmean xstd]; 
    x = x*scaling_factor;
else
    %Return all of them.
    alpha = alph_sim;
    delta = del_sim;
    lambda = lam_sim;
    x = x*scaling_factor;
end

if disp_waitbar, close(f), end

end

