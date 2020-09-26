function analyze(varargin)
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
        load(strcat(foldername,'/setup.mat'));
    case 5
        analyze_folder = 0;
        nsim           = varargin{1};
        F_sim          = varargin{2};
        alph_sim       = varargin{3};
        del_sim        = varargin{4};
        lam_sim        = varargin{5};
    otherwise
        error('Wrong number of inputs.')
end
end