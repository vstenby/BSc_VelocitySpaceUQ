%% 
clear, clc, close all

load('../data/TCV/ts4.mat');

L.L0b = L0b;
L.L1E = L1E;
L.L1p = L1p;

A = transf_matrix; b = double(S_TR_invcrime_n); xtrue = double(f_TR_coarse(:));

x = NNHGS(A, b, L, 500);




