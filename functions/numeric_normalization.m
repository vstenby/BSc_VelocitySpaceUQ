function [S2hat, W2hat, factor] = numeric_normalization(Shat, What)
%Numeric Normalization according to Birgitte's Thesis
%TODO: Regularization matrix should be normalized as well if eye is not
%used.

S2hat = Shat / norm(Shat,2);
W2hat = What / norm(What,2);

factor = norm(Shat,2)/norm(What,2);
end