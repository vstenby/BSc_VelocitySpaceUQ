function [x, Active, iter, Gnrm] = GPCG(A, b, x, MaxIts, MaxGP, MaxCG, tol)
%
%  [x, iter, Gnrm] = GPCG(A, b, x, MaxIts, MaxGP, MaxCG, tol);
%
%  Gradient-Projection Conjugate-Gradient Method for solving
%
%  x* = arg min_{x>=0} .5*x'*A*x-x'*b
%
%  Nonnegatively constrained quadratic minimization.
%
%  References:
%   [1]  Jorge J. More' and Gerardo Toraldo, "On the Solution of Large
%   Quadratic 
%        Programming Problems with Bound Constraints," SIAM Journal on Optimization, 
%        1, 1991, pp. 93-113.
%
%   Input: A  -  object defining the coefficient matrix.
%          b  -  Right hand side vector
%          x  -  initial guess 
% 
%   Optional Intputs:
%     MaxIts  -  number of iterations to perform. 
%     MaxGP   -  Max gradient projection iters per outer iteration (default 1e8).
%     MaxCG   -  Max CG iters per outer iteration (default 1e8).
%        tol  -  projected gradient tolerance for stopping.
% 
%   Output:
%          x  -  solution
%       iter  -  number of outer iterations performed
%       Gnrm  -  norm of the projected gradient at each iteration
%
%
%  J. Bardsley 12-07-2007, modifications 11-11-2010
%
  if ~isnumeric(A), Afun = A.Afun; Aparams=A.Aparams; end
  if     nargin < 4, MaxIts = 100; MaxGP = 1e8; MaxCG = 1e8; tol = 1e-5;
  elseif nargin < 5, MaxGP = 1e8; MaxCG = 1e8; tol = 1e-5;
  elseif nargin < 6, MaxCG = 1e8; tol = 1e-5;
  elseif nargin < 7, tol = 1e-5; end
  if isempty(x), x = ones(size(b)); 
  else x = max(x,0); end
  
  if isnumeric(A), Ax=A*x; else Ax=feval(Afun,x,Aparams); end 
  g = Ax-b;
  J = x(:)'*(.5*Ax(:)-b(:));  % Initial value of the cost function.
  Active = (x == 0);                        % Compute active set
  pg = -g.*((1 - Active) + Active.*(g < 0));
  pgradnorm0 = norm(pg(:));

  % Iteration history storage vectors.
  [nx,ny] = size(b);
  n = nx*ny;
  Jnrm = J;
  Gnrm = pgradnorm0;

  %  GPCG ITERATIONS:
  iter = 0;
  total_iter = 0;
  outer_flag = 0;
  gp_flag = 0;
  while outer_flag == 0;
    iter = iter + 1;
    y = x;
    xold = x;
    
    % INNER GRADIENT PROJECTION (GP) ITERATIONS:
    J_diff_max = 0;
    gp_iter = 0;
    while gp_flag == 0
      gp_iter = gp_iter + 1;
      
      % Compute Cauchy step-length parameter. Then perform linesearch.
      d = -g.*((1 - Active) + Active.*(g < 0));
      if isnumeric(A), Ad = A*d; 
      else Ad=feval(Afun,d,Aparams); end 
      init_step_param = -g(:)'*d(:) / (d(:)'*Ad(:));
      [xnew,Jnew,g,ls_paramGP,ls_flag] = linesearch(y,d,g,J,A,b,init_step_param);
      if ls_flag == 4, disp(' Line search failure. '); return; end
      
      %  Update information and check GP stopping criteria.
      Active_new = (xnew == 0);
      same_active = min( Active(:) == Active_new(:) );
      J_diff = J - Jnew;
      J_diff_max = max(J_diff,J_diff_max);
      Active = Active_new;
      y = xnew;
      J = Jnew;
      if J_diff < .25*J_diff_max | same_active | gp_iter==MaxGP, 
        gp_flag = 1;
      end

    end % while gp_flag == 0

    %  INNER CG ITERATIONS:
    %  Subspace minimization using CG to approximately compute an 
    %  unconstrained minimizer of q(delx) = 0.5*delx'*A*delx + delx'*g 
    %  where A = projected Hessian, g = projected gradient.
    xold = y; gold = g; Jold  = J; Activeold = Active;
    delx  = zeros(size(g)); resid  = -(1-Active) .* g;
    
    %  CG iterations.
    J_diff_max = 0;
    cg_iter = 0;
    cg_iter0 = 0;
    cg_flag = 0;
    while cg_flag == 0
      cg_iter = cg_iter + 1;
      d  = resid;   
      rd = resid(:)'*d(:); 
      
      %  Compute new conjugate direction p.
      if cg_iter == 1, 
          p = d; 
      else
          betak  = rd/rdlast;    
          p = d +betak*p;
      end

      %  Update delx and residual.
      if isnumeric(A), Ap = A*p; 
      else Ap = feval(Afun,p,Aparams); end 
      Ap  = (1-Active).*Ap;
      alphak  = rd/(p(:)'*Ap(:));
      delx = delx + alphak*p;   
      resid = resid - alphak*Ap; 
      rdlast = rd; 
       
      %  Check for sufficient decrease in quadratic or max iter exceeded.
      %  Note that J = J(xold + delx).
      J_diff = alphak*(p(:)'*(alphak/2 * Ap(:) - resid(:)));
      J = J - J_diff;    
      J_diff_max = max(J_diff,J_diff_max);
      if J_diff <= .1*J_diff_max | cg_iter == MaxCG
        init_step_param = 1;
        [x,J,g,ls_paramCG,ls_flag] = linesearch(xold,delx,gold,Jold,A,b,init_step_param);
        if ls_flag == 4, disp(' Line search failure. '); return; end
        stepnorm = norm(x(:)-xold(:));
        Active = (x == 0);
        pg = g.*((1 - Active) + Active.*(g < 0));
        pgradnorm = norm(pg(:));
   
        %  Check CG stopping criteria.
        Binding = Active .* (g >= 0);     %  Binding set.
        same_active = min( Active(:) == Activeold(:) );
        same_AB = min(Active(:) == Binding(:));
        if ~same_AB         %  Active ~= Binding
          cg_flag = 1;      %  Stop CG iterations.
          gp_flag = 0;      %  Start GP iterations.
        elseif ~same_active %  Active = Binding while Active ~= Activeold.
          y = x;
	      gp_flag = 1;      %  Skip GP iterations.
	      cg_flag = 1;      %  Restart CG iterations.
        elseif cg_iter < MaxCG %  Active = Binding = Activeold. Continue CG.
	      J_diff_max = 0;
          cg_iter0 = cg_iter;
          %disp('Continuing CG iterations');
        elseif cg_iter == MaxCG
          cg_flag = 1;
          gp_flag = 0;
        end
      end % if J_diff ...
    end % CG iteration
    
    total_iter = total_iter + gp_iter + cg_iter;
 
    %  Output and store numerical performance information and check
    %  stopping criteria.
    Gnrm = [Gnrm, pgradnorm];
    %fprintf('It=%d J=%6.3e |s|=%6.3e |pg|=%6.3e #gp=%d #cg=%d \n',iter, J, stepnorm, pgradnorm, gp_iter, cg_iter); 
    if pgradnorm/pgradnorm0 < tol 
      outer_flag = 1; 
      %disp('GPCG: projected gradient stopping tolerance met.');
    elseif iter==MaxIts
      outer_flag = 1;
      %disp('GPCG: maximum iterations met.');
    end
  end % while outer_flag == 0

%------------------------------------------------------------------------
%  LINE SEARCH FUNCTION
%------------------------------------------------------------------------
function [xnew,Jnew,gnew,alpha,ls_flag] = linesearch(x,p,g,J,A,b,alpha_init)
% 
%  Line search algorithm for finding approximate solution to
%    min_{alpha>0} phi(alpha),
%  where phi(alpha) = J(x + alpha*p) and J is quadratic.
%
%   Input: x  -  current estimate
%          p  -  search direction
%          g  -  gradient of J at x
%          J  -  cost function
%          A  -  coefficient matrix
%          b  -  data
%          alpha_init - initial line search parameter.
% 
%   Output: 
%       xnew  -  new estimate
%       Jnew  -  J(xnew)
%       gnew  -  grad J(xnew)
%       ls_flag - reason for exiting line search
%
  if ~isnumeric(A), Afun = A.Afun; Aparams=A.Aparams; end
  
  ls_flag = 0;
  ls_iter = 0;
  phi_0 = J;
  phi_p_0 = g(:)'*p(:);
  alpha = alpha_init;
  xnew = x + alpha*p;
  
  if min(xnew(:)>=0)  %  Unconstrained min is in feasible set.
    if isnumeric(A), Axnew = A*xnew; 
    else Axnew = feval(Afun,xnew,Aparams); end  
    gnew = Axnew-b;
    Jnew = xnew(:)'*(.5*Axnew(:)-b(:));                   
    return
  else                %  Unconstrained min is not in feasible set.
    xnew = max(xnew,0);
    if isnumeric(A), Axnew = A*xnew; 
    else Axnew = feval(Afun,xnew,Aparams); end  
    gnew = Axnew-b;
    phi_alpha = xnew(:)'*(.5*Axnew(:)-b(:));                   
  end    

  %  Find largest alpha such that x + alpha*p is in feasible set.
  ii = find(p < 0);    
  t  = -x(ii)./p(ii);
  jj = find(t > 0);
  t = t(jj);
  beta_1 = min(t(:));
  
  %  Line search iteration.
  max_iter = 10;
  mu = .1;
  ls_param1 = .1;
  ls_param2 = .5;  
  while ls_flag == 0
    ls_iter = ls_iter + 1;
    % Check sufficient decrease condition.
    if phi_alpha <= phi_0 + mu*g(:)'*(xnew(:) - x(:)) 
      Jnew = phi_alpha;
      ls_flag = 3;
      return
    end
    % Minimize the quadratic which interpolates J, g and Jnew.    
    alpha_new = -.5*phi_p_0*(alpha^2/(phi_alpha - alpha*phi_p_0 - phi_0));
    % Determine new alpha.
    m = median([ls_param1*alpha,alpha_new,ls_param2*alpha]);
    alpha = max(m,beta_1);
    % Evaluate phi(alpha).
    xnew = max(x + alpha*p,0);
    if isnumeric(A), Axnew = A*xnew; 
    else Axnew = feval(Afun,xnew,Aparams); end  
    gnew = Axnew-b;
    phi_alpha = xnew(:)'*(.5*Axnew(:)-b(:));                   
    if ls_iter >= max_iter
      disp('*** Linesearch Failure: Max Line Search Iters Exceeded ***');
      ls_flag = 4;
    end
  end % while ls_flag == 0    
  return;

