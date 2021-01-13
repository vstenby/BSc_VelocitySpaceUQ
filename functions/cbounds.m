function c = cbounds(xsamp)
% Calculates the width of the credibility bounds for 
% a matrix of samples.
%
% Usage:
%	c = cbounds(x)
% Input:
%
% x: Matrix of size (N x n),
% where N is the dimension of x, and 
% n is the number of samples. 


% npixel x nsamps matrix
qlower = quantile(xsamp,0.05/2,2);
qupper = quantile(xsamp,1-(0.05/2),2);
c = qupper - qlower;

end