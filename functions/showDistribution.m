function showDistribution(X, gridinfo, caxis_limits)
%Show the distribution in X

%Set caxis if not specified.
if nargin <= 2, caxis_limits = [min(X(:)), max(X(:))]; end

if isstruct(gridinfo)
    %Reshape X 
    if any(size(X) == 1); X = reshape(X, length(gridinfo.vperp_ax), length(gridinfo.vpara_ax)); end
    
    %Display X
    imagesc(gridinfo.vpara_ax, gridinfo.vperp_ax, X); axis xy; axis image;
    caxis(caxis_limits)
    xlabel('v_{||}', 'FontSize', 10); ylabel('v_{\perp}', 'FontSize', 10);
    set(gca, 'TickLabelInterpreter', 'latex')
    colorbar();
else
    if any(size(X) == 1); X = reshape(X, gridinfo(1), gridinfo(2)); end
    
    %Dispaly X
    imagesc(X); axis xy; axis image;
    caxis(caxis_limits);
    colorbar()
end
end