.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _CompilingCanu:

Compiling Canu
==============

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Steps
-----

#. Download a source release of canu_ to a directory like ``~/project/src``

    .. code-block:: bash

        wget https://github.com/marbl/canu/archive/v1.8.tar.gz -O canu-v1.8.tar.gz

#. Unzip that sources release and load the gcc 6 module:

    .. code-block:: bash

        tar -xzvf canu-v1.8.tar.gz
        module load gcc/6.4.0
    
#. Set the compiler to the correct values:

    .. code-block:: bash
    
        export CC=`which gcc`
        export CXX=`which g++`

#. Change to the ``canu/src`` directory and compiler

    .. code-block:: bash

        cd canu/src
        make

#. This places a binary in ``/project/brant/src/canu-1.8/Linux-amd64/bin/canu`` to which we can symlink in ``~/project/shared/bin``:

    .. code-block:: bash
    
        cd $HOME/project/shared/bin
        ln -s ../src/canu-1.8/Linux-amd64/bin/canu ./