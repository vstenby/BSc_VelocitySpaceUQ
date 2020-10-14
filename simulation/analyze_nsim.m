%% Script for analyzing the nsim numerical experiment.
clear, clc, close all

%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

problem = 'small_scale';
%problem = 'large_scale';

switch problem
    case 'small_scale'
        
        %Analyze the small scale problem.
        vparamin = -4e6;
        vparamax = 4e6;
        vparadim = 40;

        vperpmin = 1e4;
        vperpmax = 4e6;
        vperpdim = 20;

        Tpara = 2e4;
        Tperp = 2e4;
        vparadrift = 5e5;

        options.Mi = 2*1.6726e-27;

        phi = [10 20 40 70 85];
        u   = [-5e6:1e5:5e6];
        ubroadening = 1;

        [vpara, vperp, gridinfo] = construct_vgrid(vparamin,vparamax,vparadim,vperpmin,vperpmax,vperpdim);
        [x_true, xinfo] = biMaxx(vpara,vperp,Tpara,Tperp,vparadrift,options);
        x_true = x_true(:);
        
        close all
        files = dir('nsimUQ/*.mat');
        files = {files.name};

        for i=1:length(files)
            load(strcat('nsimUQ/',files{i}))
            prefix = extractBetween(files{i},1,5);
            prefix = prefix{1};
            eval(strcat('x_',prefix,'=x;'))
            eval(strcat('alph_',prefix,'=alph_sim;'))
            eval(strcat('del_',prefix,'=del_sim;'))
            eval(strcat('lam_',prefix,'=lam_sim;'))

            %Analyze each of the runs
            analyze_sim(n, x, alph_sim, del_sim, lam_sim, gridinfo);
            input('Continue? ')
            close all
        end
        
        figure
        subplot(1,3,1)
        showDistribution(abs(x_n_1e2(:,1)-x_n_1e3(:,1)),gridinfo); title('Absolute diff. 1e2 and 1e3');
        subplot(1,3,2)
        showDistribution(abs(x_n_1e3(:,1)-x_n_1e4(:,1)),gridinfo); title('Absolute diff. 1e3 and 1e4');
        subplot(1,3,3)
        showDistribution(abs(x_n_1e4(:,1)-x_n_1e5(:,1)),gridinfo); title('Absolute diff. 1e4 and 1e5');
        
        r = relerr(x_true, [x_n_1e2(:,1) x_n_1e3(:,1) x_n_1e4(:,1) x_n_1e5(:,1)]);

        figure
        semilogx([1e2 1e3 1e4 1e5],r)
        
    case 'large_scale'
        
        %Construct the default grid
        [vpara, vperp, gridinfo] = construct_vgrid();
        vparadim = length(gridinfo.vpara_ax);
        vperpdim = length(gridinfo.vperp_ax);

        %Evaluate the bi-Maxwellian on this grid with default values.
        [x_true, xinfo] = biMaxx(vpara, vperp); x_true = x_true(:);
        
        close all
        files = dir('nsimUQ2/*.mat');
        files = {files.name};

        for i=1:length(files)
            load(strcat('nsimUQ2/',files{i}))
            prefix = extractBetween(files{i},1,5);
            prefix = prefix{1};
            eval(strcat('x_',prefix,'=x;'))
            eval(strcat('alph_',prefix,'=alph_sim;'))
            eval(strcat('del_',prefix,'=del_sim;'))
            eval(strcat('lam_',prefix,'=lam_sim;'))

            %Analyze each of the runs
            analyze_sim(n, x, alph_sim, del_sim, lam_sim, gridinfo);
            input('Continue? ')
            close all
        end
        
        figure
        subplot(1,2,1)
        showDistribution(abs(x_n_1e2(:,2)-x_n_1e3(:,2)),gridinfo); title('Absolute diff. 1e2 and 1e3');
        subplot(1,2,2)
        showDistribution(abs(x_n_1e3(:,2)-x_n_1e4(:,2)),gridinfo); title('Absolute diff. 1e3 and 1e4');
        
        r = relerr(x_true, [x_n_1e2(:,2) x_n_1e3(:,2) x_n_1e4(:,1)]);

        figure
        semilogx([1e2 1e3 1e4],r)
        
end



