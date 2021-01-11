clear, clc, close all
load('ts4.mat')

%L = [L0b/max(L0b(:)) ; L1E ; L1p];
L = [L1E ; L1p];
[~,~,r] = TikhNN(transf_matrix,double(S_TR_invcrime_n),logspace(-25,-10),L,'return_relerr',true,'x_true',double(f_TR_coarse(:)))