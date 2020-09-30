function [A_norm, b_norm] = error_normalization(A,b,e)
% Normalizes with e as done in AnalyticTestCase.m

%This could be done quicker.
A_norm = zeros(size(A));
for i=1:length(b)
        A_norm(i,:) = A(i,:)/e(i);
end

b_norm = b ./ e;
end