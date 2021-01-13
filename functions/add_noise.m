function[b, e] = add_noise(b0,noise_level,background_level)
% Usage:
%   b = add_noise(b0)
%   [b,e] = add_noise(b0, noise_level, background_level)
%
% This function adds noise to the right hand side. The default noise level
% is 0.01. The noise generated noise is based on some earlier code by
% Asger Schou Jacobsen, December 2014.

if nargin == 1
    noise_level = 0.01;
    background_level = 1;
elseif nargin == 2
    background_level = 1;
end

background_errorlevel = sqrt(background_level);

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
