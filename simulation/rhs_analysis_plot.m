function rhs_analysis_plot(files)

hfig = figure('units','normalized','outerposition',[0 0 1 1]);

%Find the largest relative error
rel_max = 0;
for i=1:length(files)
   load(files{i}, 'r0th', 'r1st')
   rel_max = max([rel_max ; r0th' ; r1st']);
end

for i=1:length(files)
    
       load(files{i})
       subplot(1,length(files),i)
       semilogx(alpha_relerr,r0th, 'r-')
       hold on
       semilogx(alpha_relerr,r1st, 'b-')
       ax = gca; ax.YGrid = 'on'; ax.GridLineStyle = '-'; yticks(0:0.025:1);
       legend('0th','1st', 'AutoUpdate', 'off');
       plot(optalpha_0th, minr0,'r.','MarkerSize',15);  xlabel('\alpha','FontSize',15)
       plot(optalpha_1st, minr1,'b.','MarkerSize',15);
       ylim([0 1])
       xline(CIalpha0(1),'r-.','LineWidth',3,'alpha',0.2);
       xline(CIalpha0(2),'r-.','LineWidth',3,'alpha',0.2);
       xline(mean(alphasim0),'r-','LineWidth',3,'alpha',0.9);

       xline(CIalpha1(1),'b-.','LineWidth',3,'alpha',0.2);
       xline(CIalpha1(2),'b-.','LineWidth',3,'alpha',0.2);
       xline(mean(alphasim1),'b-','LineWidth',3,'alpha',0.9);

       switch rhs
           case 1
               title('Inverse crime','FontSize',15)
           case 2
               title('No inverse crime','FontSize',15)
           case 3
               title('Analytic right hand side','FontSize',15)
       end
       hold on
    end
end