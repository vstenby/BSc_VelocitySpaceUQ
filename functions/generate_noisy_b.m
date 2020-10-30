function [b, e] = generate_noisy_b(A, x, varargin)
% Generates a noisy right hand side.
%
% Usage: 
%    ``[b, e] = generate_noisy_b(A, x)``
%
%    ``[b, e] = generate_noisy_b(A, x, varargin)``
%
% Inputs:
%    * **A**:                System matrix
%
%    * **x**:                True solution
%
% Optional inputs:
%    * **noise_level**:      Noise level
%
%    * **background_level**: Background level for the noise.
%
% Output:
%    * **b**:                Write some description here.
%
%    * **e**:                Write some description here.

x = x(:);

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
%Set default values for optional parameters.
noise_level = 0.01;
background_level = 1;
%Unpack the varargin and evaluate.
validvars = {'noise_level','background_level'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end
% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

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