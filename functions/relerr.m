function [r, min_idx] = relerr(x_true, x)
%Finds the relative error between the true solution, x_true, and x.
%Both should be vectors.

if size(x_true,2)~=1
    x_true = reshape(x_true, numel(x_true), []);
end

nx = size(x,2);

r = vecnorm(repmat(x_true,[1,nx])-x)./norm(x_true);

[~,min_idx] = min(r);
end