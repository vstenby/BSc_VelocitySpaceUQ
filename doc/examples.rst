Examples
========

Tikhonov Reconstructions
------------------------

Bi-Maxwellian 
^^^^^^^^^^^^^

.. code-block:: Matlab

   %% Bi-Maxwellian Reconstruction Example
   % This is a demo of a velocity-tomography example, where 
   % we reconstruct the bi-Maxwellian fast-ion velocity distribution
   % with the 0th and 1st order Tikhonov formulation.

   clear, clc, close all

   [A, b, x, L, ginfo] = biMax();

   %Display the true solution
   figure
   showDistribution(x,ginfo); title('True solution')

   %Display nonnegative least-squares solution
   figure
   showDistribution(TikhNN(A,b,0,[]),ginfo); title('NNLS solution')

   %0th order Tikhonov.
   disp('Solving 0th order Tikhonov.')
   alpha0 = logspace(-15,-8,20);
   xalpha0 = TikhNN(A, b, alpha0, []);
   [ralpha0, idx0] = relerr(x, xalpha0);
   fprintf('Optimal solution: alpha = %.2e, r(alpha) = %.5f\n',alpha0(idx0),ralpha0(idx0))

   figure; 
   semilogx(alpha0, ralpha0); xlabel('\alpha'); ylabel('Relative error')
   hold on
   plot(alpha0(idx0),ralpha0(idx0), '.', 'MarkerSize',15); ylim([0 max(ralpha0)])
   title('Relative error, 0th order Tikhonov')

   figure
   showDistribution(xalpha0(:,idx0),ginfo); title('Optimal 0th order Tikhonov solution')

   %1st order Tikhonov.
   disp('Solving 1st order Tikhonov.')
   alpha1 = logspace(-4,4,20); 
   xalpha1 = TikhNN(A, b, alpha1, L);
   [ralpha1, idx1] = relerr(x, xalpha1);
   fprintf('Optimal solution: alpha = %.2e, r(alpha) = %.5f\n',alpha1(idx1),ralpha1(idx1))

   figure; 
   semilogx(alpha1, ralpha1); xlabel('\alpha'); ylabel('Relative error')
   hold on
   plot(alpha1(idx1),ralpha1(idx1), '.', 'MarkerSize',15); ylim([0 max(ralpha1)])
   title('Relative error, 1st order Tikhonov')

   figure
   showDistribution(xalpha1(:,idx1),ginfo); title('Optimal 1st order Tikhonov solution')

Slowing-down Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^

Here, we should have the slowing down distribution example.

.. code-block:: Matlab
   for i=1:5
      disp(i)
   end





Uncertainty Quantification 
--------------------------

Bi-Maxwellian Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^ 

Here is some code that shows how to perform uncertainty quantification.

Slowing-Down Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^

And here is some code that shows how to perform uncertainty quantification for the slowing-down distribution
