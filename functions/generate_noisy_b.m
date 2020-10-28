function [b, e] = generate_noisy_b(A, x, noise_level, background_level)
% [b, e] = generate_noisy_b(A,x,noise_level,background_level)
%
% This function creates a noisy b

x = x(:);

switch nargin
    case 2
        noise_level      = 0.01;
        background_level = 1;
    case 3
        background_level = 1;
end

background_errorlevel = sqrt(background_level);

%Should this be changed such that we can call biMaxb or isoSDb?
b0 = A*x;

if length(background_errorlevel) == 1
    e_min_vector = ones(size(b0))*background_errorlevel;
else
    e_min_vector = background_errorlevel;
end

err = noise_level*mean(sqrt(b0))*randn(size(b0)).*max([sqrt(b0) e_min_vector],[],2);

b = b0 + err;

% The calculated errorbar sizes calculated from the noisy signal such that
% we can not just take the square of the errorbars and cheat that way.
e = noise_level*mean(sqrt(abs(b)))*max([sqrt(abs(b)) e_min_vector],[],2);
end