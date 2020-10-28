function [x, alpha, varargout] = TikhNN(A,b,alpha,L,varargin)
% Nonnegative Tikhonov Solver.

[~,N] = size(A);

%Default parameters
solver = 'lsqnonneg'; solvers = {'GPCG','\','lsqnonneg'};
scaling = true;
nalpha = length(alpha);

%Set index for varargout.
varargout_idx = 1;

if nalpha~=1
    dispwaitbar = true;
else
    dispwaitbar = false;
end

%Make sure we deal with L.
if nargin <= 3
    L = speye(N); 
elseif isempty(L)
    L = speye(N);
end

nvarargin = length(varargin);
if mod(nvarargin,2) == 1, error('Odd number of varargin.'); end

for i=1:2:nvarargin
   switch varargin{i}
       case 'solver'
         solver = varargin{i+1};
         if ~any(strcmpi(solver,solvers)), error('Invalid solver.'), end
       case 'scaling'
         scaling = varargin{i+1};
         if ~islogical(scaling), error('Scaling should be either true or false.'), end
       case 'waitbar'
         dispwaitbar = varargin{i+1};
         if ~islogical(dispwaitbar), error('Waitbar should either be true or false.'), end
       case 'relerr'
         return_relerr = varargin{i+1};
         if ~islogical(return_relerr), error('Relative error should be either true or false.'), end
       case 'x_true'
         x_true = varargin{i+1};
       otherwise
         error('Wrong varargin.');
   end
end

%Make sure that we either have both or none
if exist('return_relerr','var') && ~exist('x_true','var')
    error("It's hard to return the relative error without the true solution")
elseif ~exist('return_relerr','var') && exist('x_true','var')
    warning('x_true given but not used.')
end


if dispwaitbar
   f = waitbar(0, 'Solving...');
end

%Scaling is set to 1/max(A) if need be, otherwise 1.
if scaling
    scaling_factor = 1/max(A(:));
else
    scaling_factor = 1;
end

%Scale A.
A = A*scaling_factor;

%Preallocate x for speed.
x = zeros(N,nalpha);

% -- Main loop --
switch solver
    case 'GPCG'
        LtL = L'*L; %Precompute LtL for speed.
        x0 = zeros(N,1); 
        for i=1:nalpha
            disp(alpha(i))
            B   = @(x) (A'*(A*x)) + alpha(i)*LtL*x; 
            rhs = A'*b;
            x(:,i) = GPCG(B, rhs, x0, 50);
            %x(:,i) = GPCG(B, rhs, x0, 50, 5, 20, 1e-6);
            %x(:,i) = GPCG(B, rhs, x0, 50, 5, 20, 1e-6); %We should be able to change this.
            if dispwaitbar, waitbar(i/nalpha, f, 'Solving...'), end
        end
    case 'lsqnonneg'
        for i=1:nalpha
            C = [A ; sqrt(alpha(i))*L];
            d = [b ; zeros(size(L,1),1)];
            x(:,i) = lsqnonneg(C,d);
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
if exist('return_relerr','var')
    x_true = reshape(x_true,numel(x_true),1); %Reshape x_true if it isn't a vector.
    varargout{varargout_idx} = relerr(x_true,x); 
    varargout_idx = varargout_idx + 1;
end

end