function [S_noisy, s] = add_noise(S0, noise_level, background_level)
%Adds noise to the samples.

if nargin == 1
    noise_level = 0.01;
    background_level = 1*10^0;
elseif nargin == 2
    background_level = 1*10^0;
elseif nargin > 3
    error('Too many inputs')
end

background_errorlevel = sqrt(background_level);
e_min_vector = ones(size(S0))*background_errorlevel;

%Elementwise standard deviation of the measurements.
s = noise_level*mean(sqrt(S0)).*max([sqrt(S0) e_min_vector],[],2); 

e = randn(size(S0)).*s;

S_noisy = S0 + e;
end