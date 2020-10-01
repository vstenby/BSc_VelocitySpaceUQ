function x = GPCG_TikhNN(A, b, alpha, L, x0)
% Solve 1/2 (||Ax-b||^2 + alpha||Lx||^2) with John Bardsley's GPCG.

[~,n] = size(A);
if nargin <= 2
    alpha = 0;
    x0 = zeros(n,1);
    L = speye(n);
elseif nargin <= 3
    x0 = zeros(n,1);
    L = speye(n);
elseif nargin <= 4
    x0 = zeros(n,1);
end

nalpha = length(alpha);
x = zeros(n,nalpha);

if nalpha > 1 
    disp_waitbar = 1;
else
    disp_waitbar = 0;
end

if disp_waitbar, f = waitbar(0, 'Starting...'); end
    
LtL = L'*L; %Precompute LtL.

for i=1:nalpha
    %Construct the B matrix.
    B = @(x) (A'*(A*x)) + alpha(i)*LtL*x; 

    %Construct the RHS, b.
    rhs = A'*b;

    %Pass this to GPCG.
    x(:,i) = GPCG(B, rhs, x0, 50, 5, 20, 1e-6); %Not sure if these should be the optimizer settings.
    
    if disp_waitbar, waitbar(i/nalpha, f, 'Solving...'), end
end
if disp_waitbar, close(f), end
end