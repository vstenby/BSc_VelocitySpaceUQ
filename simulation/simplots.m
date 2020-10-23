clear, clc, close all
%Load the neccesary functions
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

sim0th = UQSim('sim_angles_0th');
sim1st = UQSim('sim_angles_1st');

%% Initial inspection of 0th order 
figure;

i = 1;
k = 1;

while (0 < i) && (i < 18)
    subplot(1,2,1)
    showDistribution(sim0th.x(:,k,i),sim0th.gridinfo); title(num2str(sim0th.phi{i}));
    
    %Load A matrix 
    load(strcat('sim_angles_0th/',sim0th.filename{i}),'A');
    [~,S] = svd(A);
    subplot(1,2,2)
    semilogy(diag(S)); title('Singular values');
    
    k = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    switch value
        case 28
            i = i-1;
        case 29
            i = i+1;
        case 30
            k = 2;
        case 31
            k = 1;
        otherwise
    end
end

%% Initial inspection of 1st order 
figure;

i = 1;
k = 1;

while (0 < i) && (i < 18)
    subplot(1,2,1)
    showDistribution(sim1st.x(:,k,i),sim1st.gridinfo); title(num2str(sim1st.phi{i}));
    
    %Load A matrix 
    load(strcat('sim_angles_1st/',sim1st.filename{i}),'A');
    [~,S] = svd(A);
    subplot(1,2,2)
    semilogy(diag(S)); title('Singular values');
    
    k = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    switch value
        case 28
            i = i-1;
        case 29
            i = i+1;
        case 30
            k = 2;
        case 31
            k = 1;
        otherwise
    end
end

%% Simulation compare
%figure('units','normalized','outerposition',[0 0 1 1])
t = 1;
sims = {sim0th, sim1st};
simtitle = {'0th order', '1st order'};
simfiles = [17,18];
lsim = 1;
rsim = 2;
i_left  = 1;
i_right = 1;
while true
   
   %Select the sim for the left and right plot.
   simleft = sims{lsim}; simright = sims{rsim};
   
   %Grab the specific vectors.
   x_left  = simleft.x(:,t,i_left); gridinfo_left  = simleft.gridinfo;
   x_right = simright.x(:,t,i_right); gridinfo_right = simright.gridinfo;
   
   c(1) = min([x_left ; x_right]);
   c(2) = max([x_left ; x_right]);
    
   subplot(1,2,1)
   showDistribution(x_left,gridinfo_left, c); title(strcat(simtitle{lsim}, {' '}, num2str(simleft.phi{i_left})));
   subplot(1,2,2)
   showDistribution(x_right,gridinfo_right, c); title(strcat(simtitle{rsim}, {' '}, num2str(simright.phi{i_right})));
   
   k = waitforbuttonpress;
   value = double(get(gcf,'CurrentCharacter'));
   
   switch value
       
       %Toggles for the left plot!
       case 97
            %A is pressed, left plot should toggle one down.
            if i_left ~= 1
                i_left = i_left-1;
            else
                beep
            end
       case 100
            %D is pressed, left plot should toggle one up. 
            if i_left ~= simfiles(lsim)
                i_left = i_left+1;
            else
                beep
            end
           
       %Toggles for the right plot!
       case 28
            %Left arrow, right plot should toggle one down.
            if i_right ~= 1
                i_right = i_right-1;
            else
                beep
            end
       case 29
            %Right arrow, right plot should toggle one up. 
            if i_right ~= simfiles(rsim)
                i_right = i_right+1;
            else
                beep
            end
            

       case 116 %T (Toggle between mu and std)
           if t==1
               t=2;
           else
               t=1;
           end
       
       case 44 %, (Toggle simulation for the left plot)
           if lsim == 1
               lsim = 2;
           else
               lsim = 1;
           end
           
       case 46 %. (Toggle simulation for the right plot)
           if rsim == 1
               rsim = 2;
           else
               rsim = 1;
           end

       case 32  %Space bar
           close all
           break
       otherwise
           %Nothing should happen.
   end
end






















