function [x, alpha, delta, lambda, info] = NNHGS(A,b,L,n,varargin)
% Nonnegative Hierachical Gibbs Sampler

%Set the default parameters for the NNHGS
disp_waitbar = true;
welford = false;
keep = 1;
solver  = 'lsqnonneg';

%Scale our A
scaling_factor = 1/max(A(:));
A = A*scaling_factor;

nvarargin = length(varargin);

if mod(nvarargin,2) == 1
    error('Odd nuber of extra arguments')
else
    %Unpack varargin
    for i=1:2:nvarargin
       arg = varargin{i};
       switch arg
           case 'welford'
               welford = varargin{i+1};
               if ~islogical(welford)
                   error('Welford should be logical')
               end
           case 'keep'
               keep = varargin{i+1};
           case 'nburnin'
               nburnin = varargin{i+1};  
           case 'solver'
               solver = varargin{i+1};
       end
    end
end

if ~welford && keep ~= 1
   warning('Trimming has to be done manually if Welford is false') 
end

if welford
    if ~exist('nburnin','var')
        error("Burn-in should be specified with Welford's Algorithm")
    end
    %Calculate the number of 
    nidx = [1:keep:n*keep];
    nidx = nidx + nburnin;
    nsamps = nidx(end);
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

if welford
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

switch solver
    case '\'
        xtemp = [A ; sqrt(alpha0)*L]\[b ; zeros(size(L,1),1)];
    case 'lsqnonneg'
        xtemp = lsqnonneg([A ; sqrt(alpha0)*L], [b ; zeros(size(L,1),1)]);
    case 'GPCG'
        x0 = zeros(N,1);    
        B = @(x) (A'*(A*x)) + alpha0*LtL*x; 
        rhs = A'*b;
        %Should we play around with these solver settings?
        xtemp = GPCG(B, rhs, x0, 50, 5, 20, 1e-6);
    otherwise
        error('Wrong solver specified.')
end

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
if welford, welford_i = 1; end

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
   
   %bhat = mvnrnd(b,lam_temp^(-1)*speye(M))';
   %chat = mvnrnd(zeros(size(L,1),1),del_temp^(-1)*speye(size(L,1)))';
   
   switch solver
       case '\'
            %Solving (2.5) in BAHa. 
            B = [sqrt(lam_temp)*A ; sqrt(del_temp)*L];
            d = [sqrt(lam_temp)*bhat ; sqrt(del_temp)*chat];
            xtemp = B\d;
        case 'lsqnonneg'
            B = [sqrt(lam_temp)*A ; sqrt(del_temp)*L];
            d = [sqrt(lam_temp)*bhat ; sqrt(del_temp)*chat];
            xtemp = lsqnonneg(B,d);
        case 'GPCG'
            %Right hand side of (2.5) in BaHa20 divided with lambda.
            rhs = A'*bhat + alph_temp*L'*chat;
            B = @(x) (A'*(A*x)) + alph_temp*LtL*x; 
            xtemp = GPCG(B, rhs, x0, 50, 5, 20, 1e-6);
       otherwise
         error('Wrong solver specified.')
   end

   %Perhaps Welford should be written into it's own function
   if welford
       if (mod(i-nburnin,keep) || keep == 1) == 1 && nburnin < i
           xwelford = xtemp*scaling_factor; %Properly scaled of x.
           %True if i = (nburnin+1)+keep*p, where p = 0, 1, ...
           if welford_i == 1
               %M1 and S1 has to be initialized as xtemp and 0.
               Mk = xwelford;
               Sk = zeros(size(xtemp));
           else
               Mprev = Mk; %M_{k-1}
               Sprev = Sk; %S_{k-1}
               %We use the recurrence formulas
               welford_i = welford_i + 1; %so that welford_i = 2 for 2nd sample, 3 for 3rd ...
               Mk = Mprev + (xwelford - Mprev)/welford_i; 
               Sk = Sprev + (xwelford - Mprev).*(xwelford - Mk);
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

if welford
    info.welford = welford;
    info.nsamps = nsamps;
    info.nburnin = nburnin;
    info.keep = keep;
    info.nidx = nidx';
    info.alpha = alph_sim;
    info.delta = del_sim;
    info.lambda = lam_sim;
    
    %Remove the burnin and trim the three chains.
    alpha = alph_sim(info.nidx);
    delta = del_sim(info.nidx);
    lambda = lam_sim(info.nidx);
else
    %Return all of them.
    info.welford = welford;
    info.nsamps = n;
    alpha = alph_sim;
    delta = del_sim;
    lambda = lam_sim;
end

if welford
    x(:,1) = Mk; %mean
    x(:,2) = (Sk/(welford_i-1)).^(1/2); %standard deviation
else
    x = x*scaling_factor;
end

if disp_waitbar, close(f), end

end

