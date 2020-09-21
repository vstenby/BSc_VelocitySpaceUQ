function [S_noisy,e] = Generate_noisy_spectrum(transfer_matrix,distribution,noise_level,background_level)

% [S_noisy,e] = Generate_noise_spectrum(transfer_matrix,distribution,noise_level,background_level)
%
% This functions calculates a synthetic FIDA-spectrum with noise given a
% transfer_matrix, a distribution function, a certain noise_level and a
% level for the background Bremstrahlung level where the errorbars should flatten.
%
% transfer_matrix: matrix containing all the weight functions.
% distribution: column vector containing the distribution function. Should
% be in units of [fast ions pr keV pr m^3].
% noise_level: level of noise added to the signal.
% background_level: For signal levels below the background level, the
% errorbar does not scale with the signal anymore.
%
% Asger Schou Jacobsen. December 2014

if min(size(distribution)) > 1
    distribution = reshape(distribution,size(transfer_matrix,2),1);
end

background_errorlevel = sqrt(background_level);

dE = 2.5;
dp = 0.04;

%S0 = transfer_matrix*distribution*dE*dp;
S0 = transfer_matrix*distribution;

if length(background_errorlevel) == 1;
    e_min_vector = ones(size(S0))*background_errorlevel;
else
    e_min_vector = background_errorlevel;
end

S_noisy = S0 + noise_level*mean(sqrt(S0))*randn(size(S0)).*max([sqrt(S0) e_min_vector],[],2);

% Bill complained that we set the signal to 0 for negative values. He
% argues that in reality we do get negative values from noise when
% subtracting the Brehmstrahlung.

%S_noisy(S_noisy < 0) = 0;

% The calculated errorbar sizes calculated from the noisy signal such that
% we can not just take the square of the errorbars and cheat that way.
e = noise_level*mean(sqrt(abs(S_noisy)))*max([sqrt(abs(S_noisy)) e_min_vector],[],2);

figure(100); clf;
errorbar(S_noisy,e,'s')
hold all
plot(S0)
legend('Noisy signal','Noise-free signal')
% figure(110);clf;
% plot(S_noisy./e)
% title('signal-to-noise')
% figure(120);clf;
% plot(S_noisy)
% title('signal')


