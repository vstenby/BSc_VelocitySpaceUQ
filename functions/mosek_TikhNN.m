function [x, resnorm, residual, exitflag, output] = mosek_TikhNN(A, b, alphavec, L)
%Solves the regularized Tikhonov problem using mosek's lsqnonneg.
%Constructing C and d for the MOSEK routine.
%alpha should be a vector

[~,n] = size(A);
%x = zeros(n,numel(alpha));
%resnorm = zeros(numel(alpha),1);

output.algorithm = 'MOSEK';
output.iterations = zeros(length(alphavec),1);

if nargin == 3
    L = speye(n);
elseif (nargin ~= 3 && nargin~=4)
    error('Wrong number of inputs')
end

%If alpha is a vector, this should produce a waiting bar.
nalpha = length(alphavec);
if nalpha~=1
    disp_waitbar = 1;
else
    disp_waitbar = 0;
end

if disp_waitbar
    f = waitbar(0,'Solving...');
end

%TODO: Rewrite allocation of output elements to speed up the for loop.
for i=1:nalpha %If alpha is a vector.
    alpha = alphavec(i);
    
    %Mirko's original code.
    %C = [A ; sqrt(alpha) * L1vpara ; sqrt(alpha) * L1vperp ; 2e-6*sqrt(alpha)*speye(n)]; %why 2e-6?
    %d = [b ; zeros(2*size(L1vpara,2),1) ; zeros(n,1)];
    
    if isstruct(L)
        %My rewritten code to work with L1vperp and L1vpara.
        L1vperp = L.L1vperp;
        L1vpara = L.L1vpara;
        
        C = [A ; sqrt(alpha) * L1vpara ; sqrt(alpha) * L1vperp]; 
        d = [b ; zeros(2*size(L1vpara,2),1)];
    else
        %Rewritten code to work with Per Christian's L.
        C = [A ; sqrt(alpha)*L];
        d = [b ; zeros(size(L,1),1)];
    end
   
    [x(:,i), resnorm(i), residual(:,i), exitflag(i), output_temp] = lsqnonneg(C, d);
    
    if disp_waitbar, waitbar(i/nalpha, f, 'Solving...'), end
    
    if ~strcmpi(output_temp.algorithm,'MOSEK'), error('MOSEK is not used'), end
    output.iterations(i) = output_temp.iterations;
end

if disp_waitbar, close(f), end

end