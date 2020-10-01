function L = reguL(varargin)
%Make the regularization matrix L

if nargin == 2
    nvpara = varargin{1};
    nvperp = varargin{2};
else
    gridinfo = varargin{1};
    nvpara = length(gridinfo.vpara_ax);
    nvperp = length(gridinfo.vperp_ax);
end

I = @(n) speye(n);
D = @(n) spdiags([-ones(n,1) ones(n,1)],[-1 0],n+1,n);

nx = nvpara;
ny = nvperp;

D = @(n) spdiags([-ones(n,1) ones(n,1)],[-1 0],n,n);
L = [kron(D(ny), speye(nx)) ; kron(speye(ny),D(nx))];
end