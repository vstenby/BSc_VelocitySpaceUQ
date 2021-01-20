function [A, b, L, scaling] = norm_normalization(A, b, L)
%norm_normalization implements the normalization used in Birgitte's 
%Thesis. It is used for the real-data regularisation cases.

if isstruct(L)
    % Fetch from L struct
    L1p = L.L1p;
    L1E = L.L1E;
    L0b = L.L0b;
    L = [L1E ; L1p];
    norm_b = norm(b);
    b = b/norm_b;
    A = A/norm_b;
    norm_A = norm(A);
    A = A/norm_A;
    norm_L = norm(L);
    L = L/norm_L * norm(A);
    L0b = L0b/norm(L0b);
    L0b = sqrt(1e10)*L0b;
    L = [L ; L0b];
    scaling = 1/norm_A;
else
    norm_b = norm(b);
    b = b/norm_b;
    A = A/norm_b;
    norm_A = norm(A);
    A = A/norm_A;
    norm_L = norm(full(L));
    L = L/norm_L * norm(A);
    scaling = 1/norm_A;
end

end