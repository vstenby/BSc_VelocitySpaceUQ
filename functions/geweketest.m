function [p, pall] = geweketest(v, burnin, trim)
% Calculates the Geweke test statistic. Here, 0 means that the chain is
% stationary and 1 means that the chain is not.
%
% Usage: 
%    ``p = geweketest(v)``
%
%    ``[p, pall] = geweketest(v, burnin, trim)``
%
% Inputs:
%    * **v**:                Matrix where each column is a chain, e.g. ``[alpha lambda delta]``.
%
%    * **burnin**:           Removes first ``burnin`` rows, default is 0.
%
%    * **trim**:             Keep every ``trim`` elements. Default is 1.
%
% Output:
%    * **p**:                Geweke p-score for each chain.
%
%    * **pall**:             Geweke p-score for all chains.

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
    pall(i) = 1-pall(i);           %therefore we say 1-p.
end

%we take the largest p of all, i.e. the "most unstationary"-chain p-value.
p = max(pall);

end