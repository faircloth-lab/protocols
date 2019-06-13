.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _CompilingArksAndLinks:

Compiling Arks And Links
========================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Steps
-----

#. Download and unpack arks_

    .. code-block:: bash

        wget -O arks-v1.0.4.tar.gz https://github.com/bcgsc/arks/archive/v1.0.4.tar.gz
        tar -xzvf arks-v1.0.4.tar.gz

#. Download and compile a local copy of sparsehash

    .. code-block:: bash
    
        # in $HOME/project/shared/src
        cd sparsehash/ && ./configure && make && cd ../

#. Load the Boost and GCC 6 modules and set ``CC`` and ``CXX`` and ``CPPFLAGS``

    .. code-block:: bash

        module load boost/1.63.0/INTEL-18.0.0
        module load gcc/6.4.0

        export CC=`which gcc`
        export CXX=`which g++`
        export CPPFLAGS=-I$(readlink -f sparsehash/src)

#. Enter the arks_ directory and run ``autogen.sh`` and ``configure`` and, finally, ``make install``

    .. code-block:: bash
    
        cd arks-1.0.4/
        ./autogen.sh
        ./configure --with-boost=/usr/local/packages/boost/1.63.0/INTEL-18.0.0/lib --prefix=$HOME/project/shared
        make install

#. Apparently, one of the needed files does not get copies to ``bin``, so:

    .. code-block:: bash
    
        cd ~/project/shared/bin
        ln -s ../src/arks-1.0.4/Examples/makeTSVfile.py

#. Download links_ and unzip it

    .. code-block:: bash
    
        wget https://github.com/bcgsc/LINKS/releases/download/v1.8.7/links_v1-8-7.tar.gz
        tar -xvf links_v1-8-7.tar.gz

#. We need to build the bloomfilter module for links_. To compile with more modern versions of GCC (> v4), omit the ``-Dbool=char`` flag:

    .. code-block:: bash

        module load perl/5.24.0/INTEL-18.0.0
        cd links_v1.8.7/lib
        rm -rf bloomfilter/
        git clone git://github.com/bcgsc/bloomfilter.git
        cd bloomfilter/swig
        swig -Wall -c++ -perl5 BloomFilter.i
        # omit the -Dbool=char flag in the following
        g++ -c BloomFilter_wrap.cxx -I/usr/local/packages/perl/5.24.0/INTEL-18.0.0/lib/5.24.0/x86_64-linux-thread-multi/CORE -fPIC -Dbool=char -O3
        g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so -O3
        
        # test with
        perl test.pl

#. Finally, create a conda environment that contains the remaning dependencies or arks_:

    .. code-block:: bash
    
        conda create -n scaffolding tigmint bwa samtools bedtools