function L = reguL(nvpara, nvperp)
%Make the regularization matrix L

I = @(n) speye(n);
D = @(n) spdiags([-ones(n,1) ones(n,1)],[-1 0],n+1,n);

%n_x = nvpara;
%n_y = nvperp;

%D = spdiags([-ones(n_x,1) ones(n_x,1)],[-1 0],n_x+1,n_x); I = speye(n_y,n_y); 
%L1x = kron(I,D'*D);

L1vpara = kron(I(nvperp),D(nvpara)'*D(nvpara));
L1vperp = kron(D(nvperp)'*D(nvperp),I(nvpara));

%D = spdiags([-ones(n_y,1) ones(n_y,1)],[-1 0],n_y+1,n_y); I = speye(n_x,n_x); 
%L1y = kron(D'*D,I);

Lpch = [L1vpara ; L1vperp];

L = Lpch' * Lpch;


end