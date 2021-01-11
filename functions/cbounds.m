function c = cbounds(xsamp)
% npixel x nsamps matrix
qlower = quantile(xsamp,0.05/2,2);
qupper = quantile(xsamp,1-(0.05/2),2);
c = qupper - qlower;
end