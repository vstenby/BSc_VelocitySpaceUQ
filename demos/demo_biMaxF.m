% demo_biMaxF
% This script shows the different ways of using the biMaxF function
% for generating the drifting bi-Maxwellian distribution.
clear, clc, close all

fprintf('Starting demo_biMaxF:\n\n')

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%If not specified, biMaxF will have several default values. These default
%values are specified in the help command.
help biMaxF

%Also worth noting - these default values are the same as in biMaxS.

fprintf('1: Dimensions of vpara and vperp\n\n')
%biMaxF can either be evaluated on (vpara, vperp) with them being vectors
disp('Example 1')
vpara = [-1.4e7 ; 1.4e7];
vperp = [1e5 ; 1.4e7];
F = biMaxF(vpara, vperp);
fprintf('size of vpara : (%d x %d)\n', size(vpara,1), size(vpara,2))
fprintf('size of vperp : (%d x %d)\n', size(vperp,1), size(vperp,2))
fprintf('size of F     : (%d x %d)\n\n', size(F,1), size(F,2))

%or matrices (such as grids)
disp('Example 2')
[vpara, vperp] = meshgrid(linspace(-1.4e7, 1.4e7,3), linspace(1e5, 1.4e7, 3));
F = biMaxF(vpara, vperp);
fprintf('size of vpara : (%d x %d)\n', size(vpara,1), size(vpara,2))
fprintf('size of vperp : (%d x %d)\n', size(vperp,1), size(vperp,2))
fprintf('size of F     : (%d x %d)\n\n', size(F,1), size(F,2))

fprintf('2: Number of inputs and the info output\n\n')
%In terms of inputs, there are 4 options:
%1: biMaxF(vpara, vperp) 
%2: biMaxF(vpara, vperp, Tpara, Tperp)
%3: biMaxF(vpara, vperp, Tpara, Tperp, vparadrift)
%4: biMaxF(vpara, vperp, Tpara, Tperp, vparadrift, options)

disp('Example 1')
%In the options structure, we can specify Mi and ne.
options.Mi = 5*(1.6726e-27);
[F, info] = biMaxF(vpara, vperp, 200e3, 200e3, 5e6, options);
%Here, because ne is not in the options structure, it defaults to:
fprintf('Mi is specified, so Mi = %e\nne is not specified, so ne = %e\n\n', info.Mi, info.ne)

disp('Example 2')
%This can be changed if we specify it.
options.ne = 1e20;
[~, info] = biMaxF(vpara, vperp, 200e3, 200e3, 5e6, options);
fprintf('ne is now %e\n',info.ne)

%The info structure contains a lot of information about how X is generated.
%vpara      : info.vpara
%vperp      : info.vperp
%Tpara      : info.Tpara
%Tperp      : info.Tperp
%vparadrift : info.vparadrift
%Mi         : info.Mi
%ne         : info.ne