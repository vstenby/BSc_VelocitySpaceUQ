function [x, alpha, varargout] = TikhNN(A,b,alpha,L,varargin)
% Solves the Tikhonov problem formulated as:
%   
% :math:`\mathrm{x}_{\alpha} = \underset{x}{\min} \frac{1}{2} \left\| \mathrm{A} \mathrm{x}
% - \mathrm{b} \right\|^2 + \frac{\alpha}{2} \|\mathrm{L}\mathrm{x}\|^2`
% 
% or equivalently
%
% :math:`\mathrm{x}_{\alpha} = \underset{x}{\min} \frac{1}{2} \left\| \begin{bmatrix} \mathrm{A} \\ \sqrt{\alpha} \, \mathrm{L}
% \end{bmatrix} \mathrm{x} - \begin{bmatrix} \mathrm{b} \\ 0 \end{bmatrix}
% \right\|^2`
%
%
% Usage: 
%    ``[x, alpha] = TikhNN(A,b,alpha)``
%
%    ``[x, alpha] = TikhNN(A,b,alpha,L)``
%
%    ``[x, alpha, relerr] = TikhNN(A,b,alpha,L,'return_relerr',true,'xtrue',xtrue)``
%
% Inputs:
%    * **A**:               System matrix
%
%    * **b**:               Right hand side 
%
%    * **alpha**:           Regularization parameter
%
%    * **L**:               Regularization matrix
%
% Optional inputs:
%    * **solver**:          String specifying the used solver. Either ``lsqnonneg`` (default) , ``GPCG`` or ``\``.
%   
%    * **scaling**:         Whether or not A should be scaled.
%
%    * **dispwaitbar**:     Whether or not to display waitbar if several alphas are passed. 
%
%    * **return_relerr**:   If ``true``, then relerr is returned as the third output.
%
%    * **xtrue**:           Used to calculate the relative error.
%
% Output:
%    * **x**:               The Tikhonov solution for a given alpha.
%
%    * **alpha**:           Corresponding value of alpha.
%
%    * **relerr**:          Relative error

[M,N] = size(A);

%Make sure we deal with L.
if nargin <= 3
    L = speye(N); 
elseif isempty(L)
    L = speye(N);
end

%Set index for varargout.
varargout_idx = 1;

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 
%Set default values for optional parameters.
solver = 'lsqlin'; solvers = {'lsqlin','lsqnonneg','\'};
scaling = true; %Default scaling is on.
return_relerr = false;


nalpha = length(alpha);
if nalpha~=1
    dispwaitbar = true;
else
    dispwaitbar = false;
end

%Unpack the varargin and evaluate.
validvars = {'return_relerr','scaling','solver','x_true'};
evals = varargin_to_eval(varargin,validvars);
for i=1:length(evals); eval(evals{i}); end


if ~any(strcmpi(solver,solvers)), error('Invalid solver.'), end
%if ~islogical(scaling), error('Scaling should be either true or false.'), end
if ~islogical(dispwaitbar), error('Waitbar should either be true or false.'), end
if ~islogical(return_relerr), error('Relative error should be either true or false.'), end
    
%Make sure that we either have both or none
if return_relerr && ~exist('x_true','var')
    error("It's hard to return the relative error without the true solution")
elseif ~return_relerr && exist('x_true','var')
    warning('x_true given but not used.')
end

if ~check_mosek() && strcmpi(solver,'lsqnonneg')
   warning('mosek is not installed, switching default solver to lsqlin.') 
   solver = 'lsqlin';
end

if strcmpi(solver,'lsqlin')
    opts = optimset('display','off'); 
end

if isstruct(L) && ~strcmpi(scaling,'norm')
    %If L is a struct containing L1p, L1E and L0b, then 
    %we should Birgitte's scaling method.
    scaling = 'norm';
end

% - - - - - - - - - - -  Optional inputs - - - - - - - - - - - 

if dispwaitbar
   f = waitbar(0, 'Solving...');
end

%Scaling is set to 1/max(A) if need be, otherwise 1.
if strcmpi(scaling,'1/maxA')
    scaling_factor = 1/max(A(:));
    %Scale A.
    A = A*scaling_factor;
elseif strcmpi(scaling,'norm')
    [A, b, L, scaling_factor] = norm_normalization(A,b,L);
else
    scaling_factor = 1;
end



%Preallocate x for speed.
x = zeros(N,nalpha);

% -- Main loop --
switch solver
    case 'lsqnonneg'
        for i=1:nalpha
            C = [A ; sqrt(alpha(i))*L];
            d = [b ; zeros(size(L,1),1)];
            x(:,i) = lsqnonneg(C,d);
            if dispwaitbar, waitbar(i/nalpha, f, 'Solving...'), end
        end
    case 'lsqlin'
        for i=1:nalpha
            C = [A ; sqrt(alpha(i))*L];
            d = [b ; zeros(size(L,1),1)];
            lb = zeros(N,1); %Lower bound;
            [x(:,i), ~, ~, exitflag] = lsqlin(C,d,[],[],[],[],lb,[],zeros(N,1),opts);
            if exitflag <= 0, warning('not optimal solution'), end
            if dispwaitbar, waitbar(i/nalpha, f, 'Solving...'), end
        end
    case '\'
        for i=1:nalpha
            C = [A ; sqrt(alpha(i))*L];
            d = [b ; zeros(size(L,1),1)];
            x(:,i) = C\d;
            if dispwaitbar, waitbar(i/nalpha, f, 'Solving...'), end
        end
end

%Rescale back again.
x = x*scaling_factor;
if dispwaitbar, close(f), end

%Return the relative error.
if return_relerr
    x_true = reshape(x_true,numel(x_true),1); %Reshape x_true if it isn't a vector.
    varargout{varargout_idx} = relerr(x_true,x); 
    varargout_idx = varargout_idx + 1;
end
end