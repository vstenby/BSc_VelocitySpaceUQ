function [T0, T1] = phi_analysis_table(files)
%
%
%
%

% We assume that nburnin and nsamps are the same for 0th and 1st order 
% Tikhonov.

nfiles = length(files);
files = reshape(files,nfiles,[]); 
empty_cell = cell(nfiles,1);

phi_T0 = empty_cell;
geweke_p0 = empty_cell;
nburnin0 = empty_cell;
nsamples0 = empty_cell;

opt_alpha0 = empty_cell;
mean_alpha0 = empty_cell;
q_lwr0 = empty_cell;
q_upr0 = empty_cell;
qstatus0 = empty_cell;

opt_relerr0 = empty_cell;
relerr_x_posteriormean0 = empty_cell;
relerr_x_meanalpha0 = empty_cell;

phi_T1 = empty_cell;
geweke_p1 = empty_cell;
nburnin1 = empty_cell;
nsamples1 = empty_cell;

opt_alpha1 = empty_cell;
mean_alpha1 = empty_cell;
q_lwr1 = empty_cell;
q_upr1 = empty_cell;
qstatus1 = empty_cell;

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
   %CI_lwr0{i} = CIalpha0(1); 
   %CI_upr0{i} = CIalpha0(2);
   
   try
       q_lwr0{i} = qalpha0(1); 
       q_upr0{i} = qalpha0(2);
   catch
       q_lwr0{i} = CIalpha0(1); 
       q_upr0{i} = CIalpha0(2);   
   end
   
   %Fill out CI status
   if opt_alpha0{i} > q_upr0{i}
       qstatus0{i} = 'Too little';
   elseif opt_alpha0{i} < q_lwr0{i}
       qstatus0{i} = 'Too much';
   else
       qstatus0{i} = 'Good'; 
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
   
   try
       q_lwr1{i} = qalpha1(1); 
       q_upr1{i} = qalpha1(2);
   catch
       q_lwr1{i} = CIalpha1(1); 
       q_upr1{i} = CIalpha1(2);   
   end
   
   %Fill out CI status
   if opt_alpha1{i} > q_upr1{i}
       qstatus1{i} = 'Too little';
   elseif opt_alpha1{i} < q_lwr1{i}
       qstatus1{i} = 'Too much';
   else
       qstatus1{i} = 'Good'; 
   end
     
   opt_relerr1{i} = relerr(xtrue(:),xopt1);
   relerr_x_posteriormean1{i} = relerr(xtrue(:),xNNHGS1(:,1));
   relerr_x_meanalpha1{i} = relerr(xtrue(:),xsamplealpha1);
end

varnames = {'Gewkeke p-value', 'Burnin', 'Samples', ... 
            'alpha (opt)', 'alpha (mean)', '.025q', '.975q', 'Status', ...
            'Relative error (opt)', 'Relative error (sample mean)', 'Relative error (mean alpha)'};

%Round the alpha values.
opt_alpha0 = cellfun(@(x)sprintf('%.1e',x),opt_alpha0,'UniformOutput',false);
mean_alpha0 = cellfun(@(x)sprintf('%.1e',x),mean_alpha0,'UniformOutput',false);
q_lwr0 = cellfun(@(x)sprintf('%.1e',x),q_lwr0,'UniformOutput',false);
q_upr0 = cellfun(@(x)sprintf('%.1e',x),q_upr0,'UniformOutput',false);

T0 = table(geweke_p0, nburnin0, nsamples0, ...
           opt_alpha0,mean_alpha0,q_lwr0,q_upr0, qstatus0, ...
           opt_relerr0,relerr_x_posteriormean0,relerr_x_meanalpha0, ...
           'RowNames',phi_T0);
       
T0.Properties.VariableNames = varnames;

%Round the alpha values.
opt_alpha1 = cellfun(@(x)sprintf('%.1e',x),opt_alpha1,'UniformOutput',false);
mean_alpha1 = cellfun(@(x)sprintf('%.1e',x),mean_alpha1,'UniformOutput',false);
q_lwr1 = cellfun(@(x)sprintf('%.1e',x),q_lwr1,'UniformOutput',false);
q_upr1 = cellfun(@(x)sprintf('%.1e',x),q_upr1,'UniformOutput',false);

T1 = table(geweke_p1, nburnin1, nsamples1, ...
           opt_alpha1,mean_alpha1,q_lwr1,q_upr1, qstatus1, ...
           opt_relerr1,relerr_x_posteriormean1,relerr_x_meanalpha1, ...
           'RowNames',phi_T1);
       
T1.Properties.VariableNames = varnames;

clearvars -except T0 T1 

figure('Name','0th order Tikhonov Table')
uitable('Data',T0{:,:},'ColumnName',T0.Properties.VariableNames,...
    'RowName',T0.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

figure('Name','1st order Tikhonov Table')
uitable('Data',T1{:,:},'ColumnName',T1.Properties.VariableNames,...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
end