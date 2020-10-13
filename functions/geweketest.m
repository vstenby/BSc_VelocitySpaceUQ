function [p, pall] = geweketest(v)
%   Returns the p test statistic from the Geweke test.
%
%   geweketest(v)
%   where v is a N x T-dimensional matrix, where the
%   geweke test is run on each column and the largest test statistic
%   is returned.
%

T = size(v,2);
pall = zeros(T,1);

for i=1:T
    [~, pall(i)] = geweke(v(:,i)); %Returns 1 if stationary, 0 if not.
    pall(i) = 1-pall(i); 
end

p = max(pall);

end