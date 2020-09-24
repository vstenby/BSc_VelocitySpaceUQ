function [x_sim, del_sim, lam_sim, alph_sim] = NNHGS_UQ(A, b, alpha0, n)
%Nonnegative Hierachical Gibbs Sampler.

[M,N] = size(A);

%Allocate vectors for sampling.
del_sim    = zeros(n,1);
lam_sim    = zeros(n,1);
alph_sim   = zeros(n,1);
x_sim      = zeros(N,n);

%Find the initial xalpha
xtemp = mosek_TikhNN(A,b,alpha0);

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

%Calculate this beforehand.
Atb = A'*b;

for i=2:n
   disp(i)
   Axtemp = A*xtemp;
  
   %Note here that 1./ is because of the way MATLAB does gamrnd. 
   %Equation (58) and (59) in Dental CT.
   lam_sim(i) = gamrnd(a0 + M/2, 1./(t0+norm(Axtemp(:)-b(:))^2/2)); 
   del_sim(i)  = gamrnd(a1 + N/2, 1./(t1+norm(xtemp)^2/2));   
   
   alph_sim(i) = del_sim(i)/lam_sim(i);
   
   %2nd term of RHS, (69) in Dental CT.
   eta = sqrt(lam_sim(i))*A'*randn(M,1) + sqrt(del_sim(i)).*randn(N,1);
   
   %RHS is divided with lambda.
   c = Atb + eta/lam_sim(i);
   
   %LHS is divided with lambda as well such that we get alpha instead.
   B = @(x) (A' * (A*x)) + alph_sim(i)*x;
   
   xtemp = GPCG(B, c, zeros(N,1), 50, 5, 20, 1e-6);
   x_sim(:,i) = xtemp;
    
end
end