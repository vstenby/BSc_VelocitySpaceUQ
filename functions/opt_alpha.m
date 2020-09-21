function [xalpha, error_dist] = opt_alpha(alphavec, A, b, x, factor)
%Find the optimal solution given a vector of alpha values.

if nargin < 5
    factor = 1;
end

error_dist = zeros(length(alphavec),1);
xalpha = zeros(numel(x),length(alphavec));

for i=1:length(alphavec)
   xalpha(:,i) = mosek_TikhNN(A,b,alphavec(i)).*factor;
   error_dist(i) = norm(xalpha(:,i) - x(:),2);
end
end