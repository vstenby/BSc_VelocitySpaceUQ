classdef UQSim
   %Uncertainty Quantification Simulation
   
   properties
      foldername %Name of the folder with the simulations.
      x_true     %True solution
      x          %Reconstructions
      alpha      %delta-chains
      delta      %delta-chains
      lambda     %lambda-chains
      u          %u vector
      phi        %observation angles
      gridinfo   %gridinfo
      filename   %filename
   end
   
   methods 
       function obj = UQSim(foldername, varargin)
           if ~isfolder(foldername), error('Folder cannot be located.'), end
           
           if mod(length(varargin),2) == 1; error('Odd varargin.'), end;
           
           for i=1:2:length(varargin)
              switch varargin{i}
                  case 'drop'
                      drop = varargin{i+1};
                  case 'keep1st'
                      keep1stvar = varargin{i+1};
                      if ~iscell(keep1stvar); keep1stvar = {keep1stvar}; end
              end
           end
           
           if ~exist('drop','var'), drop = {''}; end
           obj.foldername = foldername;  
           [x, x_true, alpha, delta, lambda, u, phi, gridinfo, filename] = load_simulation(foldername, drop);
           obj.x_true = x_true{1};     %We assume that the true solution is the same across simulations.
           
           %Reshape x
           N = size(x{1},1);
           X = zeros(N,2,length(filename));
           
           for i=1:length(filename)
              X(:,:,i) = x{i};
           end
          
           obj.x = X;
           obj.alpha = alpha;
           obj.delta = delta;
           obj.lambda = lambda;     
           obj.gridinfo = gridinfo{1}; %We assume that all of the simulations are done on the same grid.
           obj.filename = filename;
           obj.phi = phi;
           obj.u = u;
           
           
       end
       
       function true(obj)
           showDistribution(obj.x_true,obj.gridinfo)
       end
       
       function slideshow(obj)
           %Makes a slideshow of all of the different means and
           %standard deviations
           figure('units','normalized','outerposition',[0 0 1 1])
           clim_mu  = [min(reshape(obj.x(:,1,:),1,[])) ...
                       max(reshape(obj.x(:,1,:),1,[]))];
       
           clim_std = [min(reshape(obj.x(:,2,:),1,[])) ...
                       max(reshape(obj.x(:,2,:),1,[]))];
                   
           nsims = length(obj.filename);
           idx = 1;
           while idx ~= 0 && idx ~= nsims+1
               filename_temp = strrep(obj.filename{idx},'_',' ');
               subplot(1,2,1)
               showDistribution(obj.x(:,1,idx),obj.gridinfo, clim_mu)
               title(strcat(filename_temp,' mean'))
               subplot(1,2,2)
               showDistribution(obj.x(:,2,idx),obj.gridinfo, clim_std)
               title(strcat(filename_temp,' std'))
               drawnow()
               k = waitforbuttonpress;
               value = double(get(gcf,'CurrentCharacter'));
               switch value
                   case 28 %Left arrow
                       idx = idx - 1;
                   case 29 %Right arrow
                       idx = idx + 1;
                   otherwise
                       close all
                       break
               end
           end  
           close all
       end
   end 
end