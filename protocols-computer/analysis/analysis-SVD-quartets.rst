.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running SVD Quartets
====================

:Author: Carl Oliveros, Brant C. Faircloth, Jessie Salter
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running SVD Quartets History`_ 

.. _Running SVD Quartets History: https://github.com/faircloth-lab/protocols/blob/master/protocols-computer/analysis/analysis-SVD-quartets.rst

Purpose
-------


Steps
-----

#. Install the stable channel of Docker following the instructions here: https://docs.docker.com/docker-for-mac/install/
#. Clone the `Dockerfile` repository

   .. code-block:: bash
   
       git clone git@github.com:faircloth-lab/dockerfiles.git

#. Change directory into the `ubuntu-14-wqmc` directory
#. Build the image for wQMC.  This will download all the stuff you need to run wQMC on 32-bit Ubuntu:

   .. code-block:: bash

       docker build -t faircloth/wqmc .

#. Now, run the docker image and mount a host directory to home directory of the container.  Here, we're mounting a directory we've created within the present working directory named `target`.   You can read from and write to this directory, so do your work in the container here:

   .. code-block:: bash

       docker run -i -t -v "$(pwd)"/target/:/home/generic/data faircloth/wqmc /bin/bash
   