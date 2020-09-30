function [xalpha, error_dist] = opt_alpha(alphavec, A, b, x, Lsq)
%Find the optimal solution given a vector of alpha values.

if nargin < 5
    Lsq = eye(numel(x));
end

error_dist = zeros(length(alphavec),1);
xalpha = zeros(numel(x),length(alphavec));
f = waitbar(0, 'Evaluating alphas');
for i=1:length(alphavec)
   xalpha(:,i) = mosek_TikhNN(A,b,alphavec(i),Lsq);
   error_dist(i) = norm(xalpha(:,i) - x(:),2);
   waitbar(i/length(alphavec),f,'Evaluating alphas')
end
close(f)
end