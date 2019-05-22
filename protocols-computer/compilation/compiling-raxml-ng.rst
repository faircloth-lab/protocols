.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _CompilingRaxmlng:

Compiling RAxML-NG
==================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Steps
-----

#. Download the source archive for `RAxML-NG`_ from github recursively so we get the correct dependencies.  Checkout 0.8.1-beta tag, but check the tags for `RAxML-NG`_ first to ensure there is not a newer version:

    .. code-block:: bash

        git clone --recursive https://github.com/amkozlov/raxml-ng
        git tag
        git checkout 0.8.1


.. _RAxML-NG: https://github.com/amkozlov/raxml-ng

#. Ensure that the correct modules are loaded for compilation of the source:

    .. code-block:: bash

        module load intel
        module load gcc/6.4.0
        module load impi/2018.0.128
        module load cmake

#. Set the correct ``CC`` and ``CXX`` variables so that they catch the correct IMPI versions after loading the modules:

    export CC=`which mpicc`
    export CXX=`which mpicxx`

#. Make sure you have checked out the 0.8.1 tag, create and enter a build directory, and make sure to set cmake correctly so that it will build the MPI version.  There is no INSTALL directory, so don't sweat that:

    .. code-block:: bash
        
        mkdir build && cd build
        cmake -DUSE_MPI=ON .. 

#. Assuming that is successful (you may get a warning about ``GTEST`` missing - that's fine - you just cannot run the tests):

    .. code-block:: bash

        make

#. Copy that binary into ``/project/brant/bin/``

#. Do the same thing for the non-MPI version:

    .. code-block:: bash
        
        mkdir build2 && cd build2
        cmake ..

#. Copy that binary into ``/project/brant/bin/``