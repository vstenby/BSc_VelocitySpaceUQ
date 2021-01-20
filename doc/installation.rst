

Installation
============

The MATLAB Toolbox
------------------

The VelocitySpaceUQ Toolbox can be downloaded directly from GitHub via the terminal by executing the following:

.. code-block:: bash
 
   git clone https://github.com/vstenby/BSc_VelocitySpaceUQ.git

or downloaded from the `GitHub repository. <https://github.com/vstenby/BSc_VelocitySpaceUQ/>`_

mosek
-----

The MATLAB Toolbox works with MATLAB's built-in optimizers, but it can take some time for some of the computations.
Therefore it is suggested that mosek is installed. For how to install mosek, you can find the guide 
`here. <https://docs.mosek.com/9.2/install/installation.html/>`_

Folder structure and other dependencies
---------------------------------------

The structure of the folder should be as follows:

.. code-block:: bash

  .
  ├── BSc_VelocitySpaceUQ
  │   ├── LICENSE
  │   ├── README.md
  │   ├── data
  │   ├── demos
  │   ├── doc
  │   ├── functions
  │   ├── images 
  │   ├── other
  │   ├── simulation
  │   └── testprobs
  └── aux
      ├── DrosteEffect-BrewerMap-ca40391
      ├── mosek
      └── table2latex.m


Remember to add the ``functions`` folder, ``testprobs`` folder and ``mosek`` folder to your path. Running the following code:

.. code-block:: Matlab

   %Returns 1 if mosek is properly installed, 0 otherwise.
   check_mosek()
    
   %Checks if the bi-Maxwellian testprob is added to path.
   [A, b, x, L] = biMax()

   %Checks if functions is added to path.
   TikhNN(A,b,0,[])


should check most functions. DrosteEffect-Brewermap is required for the colorscheme for the function ``UQmap``. The repository can be 
downloaded from its GitHub repository `here <https://github.com/DrosteEffect/BrewerMap>`_. ``table2latex.m`` is a MATLAB function used
for generating LaTeX-tables from MATLAB's tables. It can be found at the MATLAB File Exchange `at this link  <https://www.mathworks.com/matlabcentral/fileexchange/80386-table2latex>`_.
