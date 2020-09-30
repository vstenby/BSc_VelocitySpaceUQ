function [x, resnorm, residual, exitflag, output] = mosek_TikhNN(A, b, alphavec, L)
%Solves the regularized Tikhonov problem using mosek's lsqnonneg.
%Constructing C and d for the MOSEK routine.
%alpha should be a vector

[~,n] = size(A);
%x = zeros(n,numel(alpha));
%resnorm = zeros(numel(alpha),1);

output.algorithm = 'MOSEK';
output.iterations = zeros(length(alphavec),1);

%TODO: Rewrite allocation of output elements to speed up the for loop.
for i=1:length(alphavec) %If alpha is a vector.
    alpha = alphavec(i);
    if nargin == 3
        %0th order Tikhonov.
        C = [A ; sqrt(alpha) * speye(n)];
        d = [b ; zeros(n,1)];
    elseif nargin == 4
        %Mirko's original code.
        %C = [A ; sqrt(alpha) * L1vpara ; sqrt(alpha) * L1vperp ; 2e-6*sqrt(alpha)*speye(n)]; %why 2e-6?
        %d = [b ; zeros(2*size(L1vpara,2),1) ; zeros(n,1)];
        
        %My rewritten code to work with L1vperp and L1vpara.
        %C = [A ; sqrt(alpha) * L1vpara ; sqrt(alpha) * L1vperp]; 
        %d = [b ; zeros(2*size(L1vpara,2),1)]
        
        %Rewritten code to work with Per Christian's L.
        C = [A ; sqrt(alpha)*L];
        d = [b ; zeros(n,1)];
    else
        error('Wrong number of inputs')
    end
    [x(:,i), resnorm(i), residual(:,i), exitflag(i), output_temp] = lsqnonneg(C, d);
    if ~strcmpi(output_temp.algorithm,'MOSEK'), error('MOSEK is not used'), end
    output.iterations(i) = output_temp.iterations;
end


end