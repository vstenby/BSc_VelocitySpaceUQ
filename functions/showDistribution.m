function showDistribution(X, gridinfo)
%Show the distribution in X
imagesc(gridinfo.vpara_ax, gridinfo.vperp_ax, X); axis xy; axis image;
xlabel('v_{||}','FontSize',10)
ylabel('v_{\perp}','FontSize',10)
set(gca, 'TickLabelInterpreter', 'latex')
colorbar();
end