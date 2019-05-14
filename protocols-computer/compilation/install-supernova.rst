.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _InstallingSupernova:

Compiling Supernova
===================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

.. _supernova: https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/installation


Steps
-----

1. On whatever HPC system you are using (should be QB2 or another HPC system w/ a high RAM queue), download the code for `supernova`_.  As of writing this protocol, the current stable version is 2.1.

2. Navigate to a location on that system where you have sufficient space to unzip the software package.  10X handily provides their software with everything you need, so the unzipped package is rather large (~5 GB).

3. Unzip the sofware package:

    .. code-block:: bash

        tar -xzvf supernova-2.1.1.tar.gz

4. In ``~/.bashrc`` (or similar) or for the current session, update the ``$PATH`` to include the directory where we unpacked the supernova software:

    .. code-block:: bash

        export PATH=/home/brant/work/supernova-2.1.1:$PATH

5. Everything should be good to go, now.  You can test the software installation using a submission script like the following:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q workq
        #PBS -A <allocation>
        #PBS -l walltime=02:00:00
        #PBS -l nodes=1:ppn=20
        #PBS -V
        #PBS -N supernova_test
        #PBS -o supernova_test.out
        #PBS -e supernova_test.err


        export PATH=/home/brant/work/supernova-2.1.1:$PATH

        cd $PBS_O_WORKDIR
        supernova testrun --id=tiny

6. If the run succeeded, the ``supernova_test.out`` should contain, at the end, text that looks similar to:

    .. code-block:: bash

        Outputs:
        - Run summary:        /home/brant/work/supernova-assembly/tiny/outs/summary.csv
        - Run report:         /home/brant/work/supernova-assembly/tiny/outs/report.txt
        - Raw assembly files: /home/brant/work/supernova-assembly/tiny/outs/assembly


        Running onfinish handler...
        Waiting 6 seconds for UI to do final refresh.
        Pipestance completed successfully!

        Saving pipestance info to tiny/tiny.mri.tgz


    

