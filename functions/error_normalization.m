function [A_norm, b_norm] = error_normalization(A,b,e)
% This function normalises with the error by diving each element in b
% and each row in A with the 1/e.
% Normalizes with e as done in AnalyticTestCase.m

%This could be done quicker.
A_norm = zeros(size(A));
for i=1:length(b)
    A_norm(i,:) = A(i,:)/e(i);
end

b_norm = b ./ e;
end