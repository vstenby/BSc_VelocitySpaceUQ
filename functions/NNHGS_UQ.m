function [x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b, alpha0, L, n, options)
%Nonnegative Hierachical Gibbs Sampler.
%Perhaps make an option feature to this function.

if nargin <= 5
    options.waitbar = 1;
    options.welford = 0; 
end

if ~isstruct(options), error('Wrong options parsed to function.'), end

if options.welford
    %Number of burnin should be specified, otherwise we take floor(0.1*n)
    if ~isfield(options,'nburnin')
        warning("No burn-in specified. Welford's online algorithm might result in invalid mean and standard deviation if chains are not stationary.")
        options.nburnin = floor(0.1*n);
    end
    
    %if ~isfield(options,'welford_disp')
    %    %Display the mean and standard deviation while sampling.
    %    options.welford_disp = 0;
    %end
    
    welford_init = 0; %Welford's online algorithm has not been initialized.
end

[M,N] = size(A);

if all(size(L) == 0)
   L = speye(N); 
end

if options.waitbar
    f = waitbar(1/n,'Finding initial xalpha');
end

%Allocate vectors for sampling.
del_sim    = zeros(n,1);
lam_sim    = zeros(n,1);
alph_sim   = zeros(n,1);

if options.welford 
    x_sim = zeros(N,2); %where first column is mean, second column is standard deviation.
else
    x_sim = zeros(N,n); %where each column is a sample.
end

%Find the initial xalpha
xtemp = GPCG_TikhNN(A,b,alpha0,L);
%xtemp = mosek_TikhNN(A,b,alpha0,L);

LtL = L'*L;

%Set the initial values for alpha, x, delta and lambda.
alph_sim(1)    = alpha0;
x_sim(:,1)     = xtemp;
lam_sim(1)     = 1/norm(b(:)-A*xtemp(:))^2;
del_sim(1)     = alpha0 * lam_sim(1);

%Hyperpriors
a0         = 1; 
t0         = 0.0001; 
a1         = 1; 
t1         = 0.0001;

for i=2:n
   if options.waitbar, waitbar(i/n, f, 'Sampling...'), end
   Axtemp = A*xtemp;
   LtLxtemp = LtL*xtemp;
   
   %Note here that 1./ is because of the way MATLAB does gamrnd. 
   %Equation (58) and (59) in Dental CT
   lam_sim(i)  = gamrnd(a0 + M/2, 1./(t0+norm(Axtemp(:)-b(:))^2/2)); 
   del_sim(i)  = gamrnd(a1 + N/2, 1./(t1+xtemp(:)'*LtLxtemp(:)/2)); 
   
   alph_sim(i) = del_sim(i)/lam_sim(i);
   
   %Generate bhat and chat.
   bhat = mvnrnd(b,lam_sim(i)^(-1)*speye(M))';
   
   %NOTE: I tried to change the dimensions of chat here.
   %chat = mvnrnd(zeros(N,1),del_sim(i)^(-1)*speye(N))';
   chat = mvnrnd(zeros(size(L,1),1),del_sim(i)^(-1)*speye(size(L,1)))';
   
   %Right hand side of (2.5) in BaHa20 divided with lambda.
   rhs = A'*bhat + alph_sim(i)*L'*chat;
   
   %LHS is divided with lambda.
   B = @(x) (A' * (A*x)) + alph_sim(i)*LtL*x;
   
   xtemp = GPCG(B, rhs, zeros(N,1), 50, 5, 20, 1e-6);
   %xtemp = GPCG(B, rhs, zeros(N,1));
   
   if options.welford
       %Welford's Online Algorithm
       %First, we want to make sure we are past the nburnin.
       if options.nburnin < i
           %With Welford's online algorithm, we take the next i from
           %nburnin with mod(i,trim) == 0.
           if ~welford_init
               %M1 and S1 has to be initialized as xtemp and 0.
               Mk = xtemp; 
               Sk = zeros(size(xtemp));
               welford_init = 1;
               %welford_i = 1; %i for calculating rolling average, 
               % i from the loop can't be used if we're trimming.
           else
               Mprev = Mk; %M_{k-1}
               Sprev = Sk; %S_{k-1}
               %We use the recurrence formulas
               Mk = Mprev + (xtemp - Mprev)/i; %Rewrite this if we need to trim.
               Sk = Sprev + (xtemp - Mprev).*(xtemp - Mk);
               %wellford_i = wellford_i + 1;
           end
       end
   else
       %Save xsim.
       x_sim(:,i) = xtemp;
   end 
end

if options.welford
    %x_sim(:,1) is mean
    %x_sim(:,2) is standard deviation
    x_sim(:,1) = Mk;
    x_sim(:,2) = (Sk/(i-1)).^(1/2);
end
if options.waitbar, close(f), end
end