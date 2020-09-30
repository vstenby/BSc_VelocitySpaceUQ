function r = relerr(x_true, x)
%Finds the relative error between the true solution, x_true, and x.
%Both should be vectors.

if size(x_true,2)~=1, error('x_true should be a (N x 1)vector'), end

nx = size(x,2);

r = vecnorm(repmat(x_true,[1,nx])-x)./norm(x_true);
end