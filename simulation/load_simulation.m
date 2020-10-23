function [x, x_true, alpha, delta, lambda, phi, gridinfo, filename] = load_simulation(folder, drop)
% Input: 
%   folder
%   skip
%
% Output:
%   x
%   alpha
%   delta
%   lambda
%   phi
%   gridinfo

if nargin == 1; drop = {''}; end

%Read the files in the folder.
files = dir(strcat(folder,'/*.mat'));
files = {files.name};

%If the simulation file ends with _n.mat where n is some number,
%then we sort on on ascending n.
files = sort_files(files);

%Main for-loop loading the files.
if ~iscell(drop); drop = {drop}; end;
files = drop_files(files,drop); %Drop certain simulations
nfiles = length(files);

%Load the first file.
load(strcat(folder,'/',files{1}));

%Assume that these are the same across simulations.

gridinfoout = gridinfo;

%Allocate these before looping.
x_trueout   = cell(nfiles,1);
xout        = cell(nfiles,1);
alphaout    = cell(nfiles,1);
deltaout    = cell(nfiles,1);
lambdaout   = cell(nfiles,1);
phiout      = cell(nfiles,1); %cell since number of angles might vary.
filename    = cell(nfiles,1);
gridinfoout = cell(1,nfiles);

clear x alpha alphavec delta lambda

for i=1:length(files)
   load(strcat(folder,'/',files{i}),'x','alpha','delta','lambda','phi', 'x_true', 'gridinfo');
   x_trueout{i} = x_true;
   xout{i} = x;
   alphaout{i} = alpha;
   deltaout{i} = delta;
   lambdaout{i} = lambda;
   phiout{i} = phi;
   filename{i} = files{i};
   gridinfoout{i} = gridinfo;
end

%Rename the output variables
x        = xout;
alpha    = alphaout;
delta    = deltaout;
lambda   = lambdaout;
phi      = phiout;
x_true   = x_trueout;
gridinfo = gridinfoout;
end

function files_dropped = drop_files(files, drop)
ndrop = length(drop);
files_dropped = files;
for i=1:ndrop
   files_dropped = {files_dropped{~strcmpi(drop{i},files_dropped)}};
end
end

function n = extract_n(filename)
idx1 = regexp(filename,'_'); 
if isempty(idx1); n = []; return; end

%Select the last _ if there are several.
idx1 = idx1(end);

idx2 = regexp(filename,'.mat'); 
if isempty(idx2); n = []; return;  end

n = str2double(filename(idx1+1:idx2-1));
end

function files_sorted = sort_files(files)
%Auxil. function that sorts the cell array with filenames in it.
%They are by default sorted in ascending order by n.
nfiles = length(files);
ns = zeros(nfiles,1);
for i=1:nfiles
   n = extract_n(files{i});
   if isempty(n)
       error('n not found')
   else
       ns(i) = n;
   end
end

%Sort the files
[~, idxs] = sort(ns);
files_sorted = {files{idxs}};
end
