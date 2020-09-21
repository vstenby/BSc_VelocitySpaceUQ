function [Shat, What] = measurement_normalization(S, W, s)
%Normalize S and W according to Birgitte's thesis.

if any(size(s) == 1)
    %If s is a vector, it should be made into a diagonal matrix.
    sinv = 1./diag(s);
else
    sinv = inv(s);
end

Shat = sinv*S;
What = sinv*W;

end