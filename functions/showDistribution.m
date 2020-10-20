function showDistribution(X, gridinfo, caxis_limits)
%Show the distribution in X
if any(size(X) == 1)
    X = reshape(X, length(gridinfo.vperp_ax), length(gridinfo.vpara_ax));
end

if nargin <= 2
    caxis_limits = [min(X(:)), max(X(:))];
end
imagesc(gridinfo.vpara_ax, gridinfo.vperp_ax, X); axis xy; axis image;
caxis(caxis_limits)
xlabel('v_{||}','FontSize',10)
ylabel('v_{\perp}','FontSize',10)
set(gca, 'TickLabelInterpreter', 'latex')
colorbar();
end