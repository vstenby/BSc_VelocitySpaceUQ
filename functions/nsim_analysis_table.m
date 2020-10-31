function nsim_analysis_table(files)

nfiles = length(files);
files = reshape(files,nfiles,[]); 
empty_cell = cell(nfiles,1);

phi_T0 = empty_cell;
geweke_p0 = empty_cell;
nburnin0 = empty_cell;
nsamples0 = empty_cell;

opt_alpha0 = empty_cell;
mean_alpha0 = empty_cell;
CI_lwr0 = empty_cell;
CI_upr0 = empty_cell;
CIstatus0 = empty_cell;

opt_relerr0 = empty_cell;
relerr_x_posteriormean0 = empty_cell;
relerr_x_meanalpha0 = empty_cell;

phi_T1 = empty_cell;
geweke_p1 = empty_cell;
nburnin1 = empty_cell;
nsamples1 = empty_cell;

opt_alpha1 = empty_cell;
mean_alpha1 = empty_cell;
CI_lwr1 = empty_cell;
CI_upr1 = empty_cell;
CIstatus1 = empty_cell;

opt_relerr1 = empty_cell;
relerr_x_posteriormean1 = empty_cell;
relerr_x_meanalpha1 = empty_cell;

for i=1:length(files)
   %Load the corresponding folder.
   load(files{i});
  
   %Calculate values for 0th order table
   phi_T0{i} = num2str(phi);
   geweke_p0{i} = geweketest([alphasim0 deltasim0 lambdasim0]);
   nburnin0{i} = nburnin;
   nsamples0{i} = nsim;

   opt_alpha0{i} = optalpha_0th;
   mean_alpha0{i} = mean(alphasim0);
   CI_lwr0{i} = CIalpha0(1); 
   CI_upr0{i} = CIalpha0(2);
   
   %Fill out CI status
   if opt_alpha0{i} > CI_upr0{i}
       CIstatus0{i} = 'Underregularized';
   elseif opt_alpha0{i} < CI_lwr0{i}
       CIstatus0{i} = 'Overregularized';
   else
       CIstatus0{i} = 'Optimal alpha in CI'; 
   end
   
   opt_relerr0{i} = relerr(xtrue(:),xopt0);
   relerr_x_posteriormean0{i} = relerr(xtrue(:),xNNHGS0(:,1));
   relerr_x_meanalpha0{i} = relerr(xtrue(:),xsamplealpha0);
   
   %Calculate values for 1st order table
   phi_T1{i} = num2str(phi);
   geweke_p1{i} = geweketest([alphasim1 deltasim1 lambdasim1]);
   nburnin1{i} = nburnin; 
   nsamples1{i} = nsim;
   
   opt_alpha1{i} = optalpha_1st;
   mean_alpha1{i} = mean(alphasim1);
   CI_lwr1{i} = CIalpha1(1); 
   CI_upr1{i} = CIalpha1(2);
   
   %Fill out CI status
   if opt_alpha1{i} > CI_upr1{i}
       CIstatus1{i} = 'Underregularized';
   elseif opt_alpha1{i} < CI_lwr1{i}
       CIstatus1{i} = 'Overregularized';
   else
       CIstatus1{i} = 'Optimal alpha in CI'; 
   end
     
   opt_relerr1{i} = relerr(xtrue(:),xopt1);
   relerr_x_posteriormean1{i} = relerr(xtrue(:),xNNHGS1(:,1));
   relerr_x_meanalpha1{i} = relerr(xtrue(:),xsamplealpha1);
end

varnames = {'Observation angles', 'Gewkeke p-value', 'Burnin', 'Samples', ... 
            'alpha (opt)', 'alpha (mean)', 'Lower CI', 'Upper CI', 'Status', ...
            'Relative error (opt)', 'Relative error (posterior mean)', 'Relative error (mean alpha)'};
        
T0 = table(phi_T0, geweke_p0, nburnin0, nsamples0, ...
           opt_alpha0,mean_alpha0,CI_lwr0,CI_upr0, CIstatus0, ...
           opt_relerr0,relerr_x_posteriormean0,relerr_x_meanalpha0);
       
T0.Properties.VariableNames = varnames;
T0 = sortrows(T0,'Samples');

T1 = table(phi_T1, geweke_p1, nburnin1, nsamples1, ...
           opt_alpha1,mean_alpha1,CI_lwr1,CI_upr1, CIstatus1, ...
           opt_relerr1,relerr_x_posteriormean1,relerr_x_meanalpha1);

T1.Properties.VariableNames = varnames;
T1 = sortrows(T1,'Samples');

clearvars -except T0 T1 

figure('Name','0th order Tikhonov Table')
uitable('Data',T0{:,:},'ColumnName',T0.Properties.VariableNames,...
    'RowName',T0.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

figure('Name','1st order Tikhonov Table')
uitable('Data',T1{:,:},'ColumnName',T1.Properties.VariableNames,...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
end