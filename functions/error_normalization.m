function [A_norm, b_norm] = error_normalization(A,b,e)
% Divides each row in A by the corresponding element in e
% and divides b elementwise with e. Physics stuff.
%
% Usage: 
%    ``[A_norm, b_norm] = error_normalization(A, b, e)``
%
%
% Inputs:
%    * **xsamp**:               An N x n matrix, where each column is a sample.
%
%
% Output:
%    * **c**:                   The width of the credibility bounds.

%This could be done quicker.
A_norm = zeros(size(A));
for i=1:length(b)
    A_norm(i,:) = A(i,:)/e(i);
end

b_norm = b ./ e;
end