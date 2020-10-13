function [x, alpha, delta, lambda, info] = NNHGS(A,b,L,n,options)
% Nonnegative Hierachical Gibbs Sampler

if nargin <= 4
    %Options are not specified
    options.welford = 0;
    options.keep = 1;
    options.nburnin = floor(0.1*n);
    options.waitbar = 1;
end

if ~isfield(options,'saveallx'), options.saveallx = 0; end
if ~isfield(options,'waitbar'), options.waitbar = 1; end
if ~isfield(options,'keep'), options.keep = 1; end
if ~isfield(options,'welford'), options.welford = 0; end
if ~isfield(options,'alpha0'), options.alpha0 = 0; end
    
if ~isfield(options,'welford') && isfield(options,'keep')
   warning('NNHGS only trims if options.welford = 1') 
end
if ~isfield(options,'nburnin') && options.welford
    s1 = 'nburnin not specified, using nburnin = floor(0.1*n).';
    s2 = "Results from Welford's Online Algorithm might not be correct if chains are not stationary due to too small nburnin.";
    warning('foo:bar',strcat(s1,'\n',s2))
end

if ~isfield(options,'nburnin')
    if options.keep == 1
        options.nburnin = floor(0.1*n);
    else
        options.nburnin = floor((options.keep*n)/(1/0.1-1)); 
    end
end

%Number of samples to be generated in order to return n samples.
if options.welford
    nsamps = n*options.keep + options.nburnin;
    %Fix so that we don't do too many samps.
    nsamps = 1:options.keep:nsamps;
    nsamps = nsamps(end);
else
    nsamps = n;
end

%Allocate the delta, lambda and alpha
del_sim  = zeros(nsamps,1);
lam_sim  = zeros(nsamps,1);
alph_sim  = zeros(nsamps,1);

[M,N] = size(A);

if all(size(L) == 0)
   L = speye(N); 
end

if options.welford
    x = zeros(N,2);
else
    x = zeros(N,nsamps);
end

if options.waitbar
    f = waitbar(1/nsamps,'Finding initial xalpha');
end

%Initial sample.
alpha0 = options.alpha0;
xtemp = GPCG_TikhNN(A,b,alpha0,L);
alph_temp = alpha0;
lam_temp  = 1/norm(b(:)-A*xtemp(:))^2;
del_temp  = alpha0*lam_temp;

%Set the initial values for alpha, x, delta and lambda.
alph_sim(1)    = alph_temp;
lam_sim(1)     = lam_temp;
del_sim(1)     = del_temp;

%Calculate L'*L beforehand.
LtL = L'*L;

%Hyperpriors
a0 = 1; 
t0 = 0.0001; 
a1 = 1; 
t1 = 0.0001;

%Index for Welford's Online Algorithm
welford_i = 1;

for i=2:nsamps
   if options.waitbar, waitbar(i/nsamps, f, 'Sampling...'), end
   
   %NNHGS loop:
   Axtemp = A*xtemp;
   LtLxtemp = LtL*xtemp;
   
   %Note here that 1./ is because of the way MATLAB does gamrnd. 
   %Equation (58) and (59) in Dental CT
   
   lam_temp  = gamrnd(a0 + M/2, 1./(t0+norm(Axtemp(:)-b(:))^2/2)); 
   del_temp  = gamrnd(a1 + N/2, 1./(t1+xtemp(:)'*LtLxtemp(:)/2)); 
   
   alph_temp = del_temp/lam_temp;
   
   %Generate bhat and chat.
   bhat = mvnrnd(b,lam_temp^(-1)*speye(M))';
   chat = mvnrnd(zeros(size(L,1),1),del_temp^(-1)*speye(size(L,1)))';
   
   %Right hand side of (2.5) in BaHa20 divided with lambda.
   rhs = A'*bhat + alph_temp*L'*chat;
   
   %LHS is divided with lambda.
   B = @(x) (A' * (A*x)) + alph_temp*LtL*x;
   
   xtemp = GPCG(B, rhs, zeros(N,1), 50, 5, 20, 1e-6);
   
   %Perhaps Welford should be written into it's own function
   if options.welford
       if (mod(i-options.nburnin,options.keep) || options.keep == 1) == 1 && options.nburnin < i
           %True if i = (nburnin+1)+keep*p, where p = 0, 1, ...
           if welford_i == 1
               %M1 and S1 has to be initialized as xtemp and 0.
               Mk = xtemp; 
               Sk = zeros(size(xtemp));
           else
               Mprev = Mk; %M_{k-1}
               Sprev = Sk; %S_{k-1}
               %We use the recurrence formulas
               welford_i = welford_i + 1; %so that welford_i = 2 for 2nd sample, 3 for 3rd ...
               Mk = Mprev + (xtemp - Mprev)/welford_i; 
               Sk = Sprev + (xtemp - Mprev).*(xtemp - Mk);
           end
           welford_i = welford_i + 1;
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
info.n = n;

if options.welford
    info.nsamps = nsamps;
    info.nburnin = options.nburnin;
    info.keep = options.keep;
    info.nidx = [options.nburnin+1:options.keep:nsamps];
    info.alpha = alph_sim;
    info.delta = del_sim;
    info.lambda = lam_sim;
    
    %Remove the burnin and trim the three chains.
    alpha = alph_sim(info.nidx);
    delta = del_sim(info.nidx);
    lambda = lam_sim(info.nidx);
else
    %Return all of them.
    alpha = alph_sim;
    delta = del_sim;
    lambda = lam_sim;
end

if options.welford
    x(:,1) = Mk; %mean
    x(:,2) = (Sk/(welford_i-1)).^(1/2); %standard deviation
end

if options.waitbar, close(f), end

end

