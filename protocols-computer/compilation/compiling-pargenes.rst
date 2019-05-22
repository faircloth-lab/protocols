.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _CompilingPargenes:

Compiling Pargenes
==================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.


Steps
-----

1. Checkout the source code for `pargenes`_ from github

    .. code-block:: bash
    
        git clone --recursive https://github.com/BenoitMorel/ParGenes.git pargenes

.. _pargenes: https://github.com/BenoitMorel/ParGenes


2. Pargenes needs a couple of things to compile, including Cmake and MPI.  To get what we need on Supermike/Supermic:
   
    .. code-block:: bash
   
        module load intel/18.0.0
        module load gcc/6.4.0
        module load cmake/3.7.2/INTEL-18.0.0
        module load impi/2018.0.128


3. We also need to tell Cmake which compilers we want for it to use
   
    .. code-block:: bash
    
        export CC=`which gcc`
        export CXX=`which g++`


4. Now we should be able to compile:

    .. code-block:: bash

        cd pargenes
        ./install.sh

