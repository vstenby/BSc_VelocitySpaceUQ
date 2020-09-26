function analyze_sim(varargin)
% Analyzes the simulation in one of two ways:
%
% analyze(folder)
% or 
% analyze(nsim,F_sim,alph_sim,del_sim,lam_sim)
% 
switch nargin
    case 1
        analyze_folder = 1;
        folderName = varargin{1};
        %Fetches the nsim from foldername/setup.mat
        load(strcat(folderName,'/setup.mat'),'nsim', 'gridinfo'); 
        %should contain nsim and gridinfo.
    case 6
        analyze_folder = 0;
        nsim           = varargin{1};
        F_sim          = varargin{2};
        alph_sim       = varargin{3};
        del_sim        = varargin{4};
        lam_sim        = varargin{5};
        gridinfo       = varargin{6};
    otherwise
        error('Wrong number of inputs.')
end

if analyze_folder
    load_sims = 1;
    sim_number = 1; %Number of sim in folder, sim1.mat, sim2.mat ... 
    while load_sims
        try
            load(strcat(folderName,'/sim',num2str(sim_number),'.mat'))
            fprintf('Analysis of sim%d\n',sim_number)
            
            %Call the analysis function
            analyze_single_sim(nsim, F_sim, alph_sim, lam_sim, del_sim, gridinfo)
            disp(' ')
            
            val = input('Press enter to continue to next simulation or type 0 to stop. ');
            if val == 0
               load_sims = 0;
            else
                close all
                disp(' ')
                sim_number = sim_number + 1;
            end  
        catch
            load_sims = 0;
        end
    end     
else %IF NO FOLDER
    analyze_single_sim(nsim, F_sim, alph_sim, lam_sim, del_sim, gridinfo)
end
%END THE FUNCTION HERE.
end

function analyze_single_sim(nsim, F_sim, alph_sim, lam_sim, del_sim, gridinfo)
% 
% Analyze the simulation
%

nburnin = floor(0.1*nsim);

%Remove burnin
alph_sim = alph_sim(nburnin+1:end);
del_sim = del_sim(nburnin+1:end);
lam_sim = lam_sim(nburnin+1:end);
F_sim = F_sim(:, nburnin+1:end);

%Geweke tests.
%1 = stationary, 0 = not stationary.
[~, palpha] = geweke(alph_sim); 
[~, pdelta] = geweke(del_sim);
[~, plambda] = geweke(lam_sim);

%Finding our largest value in case one of the chains is NOT
%converged.
p = 1-max([palpha pdelta plambda]); 

fprintf('n = %d, nburnin = %d, geweke p = %f\n',nsim,nburnin, p)

%Calculate and display quantiles.
disp('Quantiles:')
q_alpha = print_quantiles(alph_sim,'alpha');
q_delta = print_quantiles(del_sim,'delta');
q_lambda = print_quantiles(lam_sim,'lambda');
disp(' ')

%Plot the chains
input('Press enter to show plot of alpha, lambda and delta chains ')
close all
figure(1)
plot(alph_sim); title('\alpha chain')
figure(2)
plot(del_sim); title('\delta chain')
figure(3)
plot(lam_sim); title('\lambda chain')

%Plot the histograms
input('Press enter to show histograms and scatterplot ')
close all

nbins = 50;
figure(1)
subplot(1,3,1)
hist(alph_sim,nbins); title('\alpha'); axis square
subplot(1,3,2)
hist(del_sim,nbins); title('\delta'); axis square
subplot(1,3,3)
hist(lam_sim,nbins); title('\lambda'); axis square
%histogram works weird with lambda, using hist instead.

figure(2)
scatter(del_sim,lam_sim,'.k')
xlabel('\delta')
ylabel('\lambda')


input('Press enter to show mean and std of F ')
close all
F_mu = mean(F_sim,2);
F_std = std(F_sim,0,2);

figure(1)
showDistribution(F_mu,gridinfo); title('Mean F');
figure(2)
showDistribution(F_std,gridinfo); title('Standard deviation F');
end

function q = print_quantiles(v,str)
    q = quantile(v, [0 0.25 0.5 0.75 1]);
    fprintf('%s %0.2e %0.2e %0.2e %0.2e %0.2e\n', ...
            pad(str,10), q(1), q(2), q(3), q(4), q(5))
end



