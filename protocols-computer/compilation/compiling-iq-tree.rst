.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _CompilingIQtree:

Compiling IQ-tree
=================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Steps
-----

#. Download the source archive for `IQ-Tree`_:

    .. code-block:: bash

        cd /project/brant/shared/src
        wget https://github.com/Cibiv/IQ-TREE/archive/v1.6.10.tar.gz -O iq-tree-v1.6.10.tar.gz
        unzip iq-tree-v1.6.10.tar.gz
        # this makes IQ-TREE-1.6.10 in src

.. _IQ-Tree: https://github.com/Cibiv/IQ-TREE/


#. Before compiling `IQ-Tree`_, you also need to download and compile/install `Eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ .  That's not particularly complex, except that I've found you really do need to use cmake to "install" Eigen3 to an include directory.  What follows are the steps to do that:

    .. code-block:: bash

        # put Eigen3 source into a tmp dir
        mkdir /project/brant/shared/tmp && cd /project/brant/shared/tmp

        # download Eigen3 source & unzip
        http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz -O eigen-v3.3.7.tar.gz
        tar -xzvf eigen-v3.3.7.tar.gz

        # build
        cd eigen-eigen-323c052e1731/
        mkdir build && cd build
        cmake .. -DCMAKE_INSTALL_PREFIX=/project/brant/shared/
        make install

        # this will install Eigen3 to /project/brant/shared/include/eigen3/

Difference Levels of Parallelism
::::::::::::::::::::::::::::::::

IQ-Tree, like RAxML, can use different levels of parallelism.  To achieve all the options, we need to compile each version.  There are essentially 3 choices: MPI, OMP, MPI-OMP Hybrid

OMP (only) Version
..................

#. Ensure that the correct modules are loaded for compilation of the source:

    .. code-block:: bash

        module load gcc/6.4.0
        module load cmake

#. Set the correct ``CC`` and ``CXX`` variables so that they catch the correct Intel Compiler versions after loading the modules.  Here, we're using GCC because there is currently an error using ICC on @supermike (probably because C++ library for ICC is old).

    .. code-block:: bash

        export CC=`which gcc`
        export CXX=`which g++`

#. Now, use cmake to create the makefiles and compile with make:

    .. code-block:: bash
        
        cd /project/brant/shared/src/IQ-TREE-1.6.10
        mkdir build2 && cd build2
        cmake -DEIGEN3_INCLUDE_DIR=/project/brant/shared/include/eigen3 -DIQTREE_FLAGS=omp ..
        make


MPI (only) Version
..................

#. Ensure that the correct modules are loaded for compilation of the source:

    .. code-block:: bash

        module load intel
        module load impi/2018.0.128
        module load cmake

#. Set the correct ``CC`` and ``CXX`` variables so that they catch the correct IMPI versions after loading the modules:

    .. code-block:: bash

        export CC=`which mpicc`
        export CXX=`which mpicxx`

#. Now, use cmake to create the makefiles and compile with make:

    .. code-block:: bash
        
        cd /project/brant/shared/src/IQ-TREE-1.6.10
        mkdir build && cd build
        cmake -DEIGEN3_INCLUDE_DIR=/project/brant/shared/include/eigen3 -DIQTREE_FLAGS=mpi ..
        make


MPI & OMP Hybrid Version
........................

#. Ensure that the correct modules are loaded for compilation of the source:

    .. code-block:: bash

        module load intel
        module load impi/2018.0.128
        module load cmake

#. Set the correct ``CC`` and ``CXX`` variables so that they catch the correct IMPI versions after loading the modules:

    .. code-block:: bash

        export CC=`which mpicc`
        export CXX=`which mpicxx`

#. Now, use cmake to create the makefiles and compile with make:

    .. code-block:: bash
        
        cd /project/brant/shared/src/IQ-TREE-1.6.10
        mkdir build3 && cd build3
        cmake -DEIGEN3_INCLUDE_DIR=/project/brant/shared/include/eigen3 -DIQTREE_FLAGS=omp-mpi ..
        make

