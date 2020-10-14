function [p, pall] = geweketest(v, burnin, trim)
%   Returns the p test statistic from the Geweke test.
%
%   geweketest(v)
%   where v is a N x T-dimensional matrix, where the
%   geweke test is run on each column and the largest test statistic
%   is returned.
%

if nargin == 1
   burnin = 0;
end

if nargin == 2
   %Remove burnin from v
   v = v(burnin+1:end,:);
end

if nargin == 3
   v = v(burnin+1:end,:);
   v = v(1:trim:end, :); 
end

T = size(v,2);
pall = zeros(T,1);

for i=1:T
    [~, pall(i)] = geweke(v(:,i)); %Returns 1 if stationary, 0 if not.
    pall(i) = 1-pall(i); 
end

p = max(pall);

end