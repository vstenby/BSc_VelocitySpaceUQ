function S_noisy = add_noise(S0)
%Adds noise to the samples.
noise_level=0.01;
background_level=1*10^0;
background_errorlevel = sqrt(background_level);
e_min_vector = ones(size(S0))*background_errorlevel;
S_noisy = S0 + noise_level*mean(sqrt(S0))*randn(size(S0)).*max([sqrt(S0) e_min_vector],[],2);
end