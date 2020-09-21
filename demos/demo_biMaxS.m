% demo_biMaxS
% This script shows the different ways of using the biMaxS function
% for generating the projections of the drifting bi-Maxwellian distribution.

clear, clc, close all

fprintf('Starting demo_biMaxX:\n\n')

%Add dependencies.
addpath('../functions')
addpath(genpath('../../aux'))

%If not specified, biMaxS will have several default values. These default
%values are specified in the help command.
help biMaxS

%Also worth noting - these default values are the same as in biMaxX.

fprintf('1: Ways to input u\n\n')
%There are several ways of inputting u, all constructed with inspiration
%from Mirko.

%Specify observation angles (phi)
phi=[10 20 40 70 85];

%We can simply pass u as a vector
disp('Example 1: u as vector')
u = [-3 -2 -1 0 1 2 3];
[S, info] = biMaxS(u,phi);
disp(info.u)

%We can pass u as a structure in one of the following ways:
%   u.struct containing:
%    (a) umax, udim
%    (b) umax, du
%    (c) umin, umax, udim (udim has to be an even number)
%    (d) umin, umax, du
%
%If we pass u like this, then vector u will be constructed as:
%    (a) linspace(-umax, umax, udim)
%    (b) [-umax:du:umax]
%    (c) [linspace(-umax, -umin, udim/2) linspace(umin, umax, udim/2)]
%    (d) [-umax:du:-umin umin:du:umax]

%Example of (a)
disp('Example 2: umax and udim specified')
ustruct.umax = 3; ustruct.udim = 7;
[~, info] = biMaxS(ustruct,phi);
disp(info.u)
clear ustruct

%Example of (b)
disp('Example 3: umax and du specified')
ustruct.umax = 3; ustruct.du = 1;
[~, info] = biMaxS(ustruct,phi);
disp(info.u)
clear ustruct

%Example of (c)
disp('Example 4: umax, umin and udim specified')
ustruct.umax = 3; ustruct.umin = 1; ustruct.udim = 6; %udim has to be an even number!
[~,info] = biMaxS(ustruct,phi);
disp(info.u)
clear ustruct

%Example of (d)
disp('Example 5: umax, umin and du specified')
ustruct.umax = 3; ustruct.umin = 1; ustruct.du = 1;
[~, info] = biMaxS(ustruct,phi);
disp(info.u)
clear ustruct

%To make it even more fun, instead of giving umin/umax, we can pass Emax 
%instead of umax, and Emin instead of umin.

%Example similar to (b), but with Emax instead of umax.
disp('Example 6: Emax and udim specified')
ustruct.Emax = 4e6; ustruct.udim = 7;
[~, info] = biMaxS(ustruct,phi);
disp(info.u)
clear ustruct

%Example similar to (c), but with Emax and Emin.
disp('Example 7: Emax, Emin and udim specified')
ustruct.Emax = 4e6; ustruct.Emin = 4e6; ustruct.udim = 6;
[~, info] = biMaxS(ustruct,phi);
disp(info.u)
clear ustruct

fprintf('2: Number of inputs and the info output\n\n')
%In terms of inputs, there are 4 options:
%1: biMaxS(u, phi)
%2: biMaxS(u, phi, Tpara, Tperp)
%3: biMaxS(u, phi, Tpara, Tperp, vparadrift)
%4: biMaxS(u, phi, Tpara, Tperp, vparadrift, options)

disp('Example 8: Changing values in options')
%In the options structure, we can specify Mi and ne.
options.Mi = 5*(1.6726e-27);
[~, info] = biMaxS(u, phi, 200e3, 200e3, 5e6, options);
%Here, because ne is not in the options structure, it defaults to:
fprintf('Mi is specified, so Mi = %e\nne is not specified, so ne = %e\n\n', info.Mi, info.ne)

%This can be changed if we specify it.
options.ne = 1e20;
disp('ne is set in options.')
[~, info] = biMaxS(u, phi, 200e3, 200e3, 5e6, options);
fprintf('ne is now %e\n',info.ne)

%The info structure contains a lot of information about how X is generated.
%u (vector) : info.u
%phi        : info.vperp
%Tpara      : info.Tpara
%Tperp      : info.Tperp
%vparadrift : info.vparadrift
%Mi         : info.Mi
%ne         : info.ne


