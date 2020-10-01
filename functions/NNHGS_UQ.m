function [x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b, alpha0, L, n, disp_waitbar)
%Nonnegative Hierachical Gibbs Sampler.

if nargin <= 5
    disp_waitbar = 0;  
end

[M,N] = size(A);

if disp_waitbar
    f = waitbar(1/n,'Finding initial xalpha');
end

%Allocate vectors for sampling.
del_sim    = zeros(n,1);
lam_sim    = zeros(n,1);
alph_sim   = zeros(n,1);
x_sim      = zeros(N,n);

%Find the initial xalpha
xtemp = mosek_TikhNN(A,b,alpha0,L);

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
   if disp_waitbar, waitbar(i/n, f, 'Sampling...'), end
   Axtemp = A*xtemp;
   LtLxtemp = LtL*xtemp;
   
   %Note here that 1./ is because of the way MATLAB does gamrnd. 
   %Equation (58) and (59) in Dental CT
   lam_sim(i)  = gamrnd(a0 + M/2, 1./(t0+norm(Axtemp(:)-b(:))^2/2)); 
   del_sim(i)  = gamrnd(a1 + N/2, 1./(t1+xtemp(:)'*LtLxtemp(:)/2)); 
   
   alph_sim(i) = del_sim(i)/lam_sim(i);
   
   %Generate bhat and chat.
   bhat = mvnrnd(b,lam_sim(i)^(-1)*speye(M))';
   chat = mvnrnd(zeros(N,1),del_sim(i)^(-1)*speye(N))';
   
   %Right hand side of (2.5) in BaHa20 divided with lambda.
   rhs = A'*bhat + alph_sim(i)*L*chat;
   
   %LHS is divided with lambda.
   B = @(x) (A' * (A*x)) + alph_sim(i)*LtL*x;
   
   xtemp = GPCG(B, rhs, zeros(N,1), 50, 5, 20, 1e-6);
   %xtemp = GPCG(B, rhs, zeros(N,1));
   
   x_sim(:,i) = xtemp;
    
end
if disp_waitbar, close(f), end
end