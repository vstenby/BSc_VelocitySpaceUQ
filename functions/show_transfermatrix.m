function show_transfermatrix(A, ginfo)
timedelay = 0.05;
[M, ~] = size(A);
h = figure;
caxis = [min(A(:)), max(A(:))];
for i=1:M
   showDistribution(A(i,:),ginfo, caxis); 
   title(sprintf('Row %d of %d',i,M));
   drawnow()
   pause(timedelay);
end
end