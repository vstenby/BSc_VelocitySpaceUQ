function [x, resnorm, residual, exitflag, output] = mosek_TikhNN(A, b, alpha)
%Solves the regularized Tikhonov problem using mosek's lsqnonneg.
%Constructing C and d for the MOSEK routine.
[m,n] = size(A);

C = [A ; sqrt(alpha) * eye(n)];
d = [b ; zeros(n,1)];

[x, resnorm, residual, exitflag, output] = lsqnonneg(C, d);
end