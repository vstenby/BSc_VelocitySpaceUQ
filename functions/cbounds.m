function c = cbounds(xsamp)
% Calculates the width of 95% credibility bounds. 
%
% Usage: 
%    ``c = cbounds(xsamp)``
%
%
% Inputs:
%    * **xsamp**:               An N x n matrix, where each column is a sample.
%
%
% Output:
%    * **c**:                   The width of the credibility bounds.


% npixel x nsamps matrix
qlower = quantile(xsamp,0.05/2,2);
qupper = quantile(xsamp,1-(0.05/2),2);
c = qupper - qlower;
end