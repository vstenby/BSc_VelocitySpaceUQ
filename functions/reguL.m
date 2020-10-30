function L = reguL(varargin)
%
% 
% 
switch nargin
    case 1
        dim = varargin{1}; 
        nrow = dim(1); ncol = dim(2);
    case 2
        nrow = varargin{1}; ncol = varargin{2};
    otherwise
        error('Wrong number of inputs')
end


%Constructs a regularization matrix L for a matrix X of dimensions m, n

I = @(d) speye(d);
D = @(d) spdiags([zeros(d-1,1) -ones(d-1,1) ones(d-1,1)],-1:1,d-1,d);

L = [kron(D(ncol), I(nrow)) ; kron(I(ncol),D(nrow))];

end