function [S_norm, A_norm, factor] = normalize_problem(S, A, signal_noise)
%Normalize the problem according to Birgitte Madsen's Thesis.
%Chapter 4.5: Normalizing the problem

[M, ~] = size(A);

%Assuming the signal noise is the same everywhere?
sigma_s_inv = 1/signal_noise .* eye(M);

Shat = sigma_s_inv * S;
What = sigma_s_inv * A;

S2hat = Shat ./ norm(Shat,2);
W2hat = What ./ norm(What,2);

factor = norm(Shat,2)/norm(What,2);

S_norm = S2hat;
A_norm = W2hat;
end



