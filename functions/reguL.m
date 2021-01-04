function L = reguL(varargin)
% Constructs the 1st order Tikhonov Regularization matrix without any
% assumptions about the border.
%
% Usage: 
%    ``L = reguL(vpara, vperp)`` 
%
% Inputs:
%    * **A**:                   TODO: Fix this documentation.
%
% Output:
%    * **L**:               	The regularization matrix L.


switch nargin
    case 1
        ginfo = varargin{1};
        %ginfo struct.
        s = ginfo.vperp_ax; hs = s(2)-s(1); m = length(s);
        t = ginfo.vpara_ax; ht = t(2)-t(1); n = length(t);
    case 2
        %Two grids. 
        T = varargin{1}; t = T(1,:); ht = t(2)-t(1); n = length(t);
        S = varargin{2}; s = S(:,1); hs = s(2)-s(1); m = length(s);
        
        %Check that grids are given properly
        if (ht == 0||hs == 0), error('hs or ht was zero.'), end
    otherwise
        error('Wrong number of inputs')
end

%Constructs a regularization matrix L for a matrix X of dimensions m, n
I = @(d) speye(d);
D = @(d) spdiags([zeros(d-1,1) -ones(d-1,1) ones(d-1,1)],-1:1,d-1,d);

L = [kron(1/ht*D(n), I(m)) ; kron(I(n),1/hs*D(m))];

end