clear, clc, close all

%Download the files DTUs servers.
if ~isfolder('biMax_sim_phi3')
    system('scp -r s174483@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/biMax_sim_phi3 biMax_sim_phi3');
end

%Add dependencies to path.
addpath(genpath('../functions'))
addpath(genpath('../../aux'))

addpath('biMax_sim_phi3')

files = {dir('biMax_sim_phi3/*.mat').name};

for i=1:length(files)
   file = files{i}; load(file)
   
   figure
   subplot(4,3,[1,4,7,10])
   semilogx(alpha_relerr,r0th, 'r-')
   hold on
   semilogx(alpha_relerr,r1st, 'b-')
   ax = gca; ax.YGrid = 'on'; ax.GridLineStyle = '-'; yticks(0:0.025:1);
   legend('0th','1st', 'AutoUpdate', 'off');
   plot(optalpha_0th, minr0,'r.','MarkerSize',15);  xlabel('\alpha','FontSize',15)
   plot(optalpha_1st, minr1,'b.','MarkerSize',15);
   
   xline(CIalpha0(1),'r-.','LineWidth',3,'alpha',0.2);
   xline(CIalpha0(2),'r-.','LineWidth',3,'alpha',0.2);
   xline(mean(alphasim0),'r-','LineWidth',3,'alpha',0.9);

   xline(CIalpha1(1),'b-.','LineWidth',3,'alpha',0.2);
   xline(CIalpha1(2),'b-.','LineWidth',3,'alpha',0.2);
   xline(mean(alphasim1),'b-','LineWidth',3,'alpha',0.9);

   legend('show')
   subplot(4,3,2)
   showDistribution(xopt0,ginfo,caxis_mu); title('Optimal alpha (0th)')

   subplot(4,3,3)
   showDistribution(xopt1,ginfo,caxis_mu); title('Optimal alpha (1st)')

   subplot(4,3,5)
   showDistribution(xNNHGS0(:,1),ginfo,caxis_mu); title('Posterior mean (0th order)')
   subplot(4,3,8)
   showDistribution(xsamplealpha0,ginfo,caxis_mu); title('Reconstruction with mean alpha (0th order)')
   subplot(4,3,11)
   showDistribution(xNNHGS0(:,2),ginfo,caxis_std); title('Standard deviation (0th order)')

   subplot(4,3,6)
   showDistribution(xNNHGS1(:,1),ginfo,caxis_mu); title('Posterior mean (1st order)')
   subplot(4,3,9)
   showDistribution(xsamplealpha1,ginfo,caxis_mu); title('Reconstruction with mean alpha (1st order)')
   subplot(4,3,12)
   showDistribution(xNNHGS1(:,2),ginfo,caxis_std); title('Standard deviation (1st order)')
   
   input('Press enter to continue ')
end
  