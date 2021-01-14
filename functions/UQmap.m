function UQmap(xcompare, xbaseline, ginfo, varargin)

%Set default parameters
cbound = true;
display = 'showDistribution';
%Unpack the varargin and evaluate.
validvars = {'cbound','display'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end
if ~islogical(cbound), error('cbound should be logical.'), end

if cbound
    c_compare  = cbounds(xcompare);
    c_baseline = cbounds(xbaseline);
    dif = c_compare-c_baseline; 
else
    %If Welford
    xstd_compare  = xcompare(:,2); 
    xstd_baseline = xbaseline(:,2); 
    dif = xstd_compare-xstd_baseline;
end

if strcmpi(display, 'showDistribution')
    %Some function to display x here.
    showDistribution(dif, ginfo); axis image; colorbar();
elseif strcmpi(display, 'imagesc')
    %Show it using imagesc
    imagesc(reshape(dif,ginfo(1),ginfo(2))); axis image; colorbar()
else
    error('Wrong display input')
end
caxis([-max(abs(dif)), max(abs(dif))])
try
    colormap(brewermap([],'*RdBu'))
catch
    error('brewermap not installed.')
end
end