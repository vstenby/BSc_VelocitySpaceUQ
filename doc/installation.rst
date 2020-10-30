

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

Folder structure
----------------

The structure of the folder should be as follows:

.. code-block:: bash
   
  .
  ├── BSc_VelocitySpaceUQ
  │   ├── LICENSE
  │   ├── README.md
  │   ├── demos
  │   ├── doc
  │   ├── functions
  │   ├── images
  │   ├── other
  │   ├── simulation
  │   └── testprobs
  └── aux
      └── mosek

Remember to add the functions and mosek folder to your MATLAB path before running the code.
Certain functions will give warnings if mosek is not loaded.
