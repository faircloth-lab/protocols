.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running Pargenes
================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running Pargenes`_ 

.. _Running Pargenes: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-pargenes.rst

Purpose
-------

Pargenes_ is a pipeline to generate gene trees from a large set of loci, using the most appropriate site rate substitution model.


.. _pargenes: https://github.com/BenoitMorel/ParGenes

Preliminary Steps
-----------------

#. To compile Pargenes, see :ref:`CompilingPargenes`


Steps
-----

#. Before running Pargenes_, you need to prepared your data by following several steps.  The easiest thing to do is to take the directory of loci that you wish to analyze (say, from a 75% matrix... or all loci [then subset]), and convert those loci to FASTA format:

    .. code-block:: bash

        phyluce_align_convert_one_align_to_another \
            --alignments input-alignments \
            --output input-alignments-fasta \
            --cores 12 \
            --log-path ./ \
            --input-format nexus \
            --output-format fasta


#. After formatting loci in FASTA format, you probably want to go ahead and reduce those loci, if needed, so that identical sequences for different taxa are removed.  This requires a recent version of phyluce_ (which is not, yet, publicly available).  Reduce the FASTA alignments by:

    .. code-block:: bash

        phyluce_align_reduce_alignments_with_raxml \
            --alignments input-alignments-fasta \
            --output input-alignments-fasta-reduced \
            --input-format fasta \
            --cores 12

#. After reducing your loci, you want to upload those to HPC.  Before uploading, it's probably best to package them up as ``.tar.gz`` and unpack them on Supermike/Supermic:

    .. code-block:: bash

        tar -czf input-alignments-fasta-reduced.tar.gz input-alignments-fasta-reduced
    
#. After uploading to Supermike/Supermic using something like ``rsync``, in your working directory, create a job submission script that looks like the following (be sure to use your ``<allocation>``).  This will run a "test-run" of pargenes and estimate the number of cores that we should use for optimal run-times:


    .. code-block:: bash

        #PBS -A hpc_allbirds02
        #PBS -l nodes=1:ppn=16
        #PBS -l walltime=2:00:00
        #PBS -q checkpt
        #PBS -N pargenes

        cd $PBS_O_WORKDIR
        CORES=16

        python /project/brant/shared/src/pargenes-mpi/pargenes.py \
            -a input-alignments-fasta-reduced \
            -o input-alignments-fasta-reduced-pargenes-dry-run \
            -d nt \
            -m \
            -c $CORES \
            --dry-run

#. This will produce an output folder (``input-alignments-fasta-reduced-pargenes-dry-run``).  In that folder is a log file that will contain an estimate of the number of cores we need to run a job optimally.  Remember that number.  Based on that number, setup a new ``qsub`` file for the "real" run of Pargenes_ where you adjust ``nodes=XX:ppn=YY`` and ``CORES``.  That will look something like the following, which will used 512 CPU cores to: (1) estimate the best site-rate substitution model for each locus, (2) estimate the best ML gene tree for each locus based on the most appropriate model, and (3) generate 200 boostrap replicates for everything:

    .. code-block:: bash

        #PBS -A <allocation>
        #PBS -l nodes=32:ppn=16
        #PBS -l walltime=12:00:00
        #PBS -q checkpt
        #PBS -N pargenes_dry_run

        cd $PBS_O_WORKDIR
        CORES=512

        python /project/brant/shared/src/pargenes-mpi/pargenes/pargenes-hpc.py \
            -a input-alignments-fasta-reduced \
            -o input-alignments-fasta-reduced-pargenes-bootstraps \
            -d nt \
            -m \
            -c $CORES \
            --bs-trees 200

#. Before downloading, you probably want to zip everything up, which you can do by creating a packaging ``qsub`` script like:

    .. code-block:: bash

        #PBS -A <allocation>
        #PBS -l nodes=1:ppn=16
        #PBS -l walltime=6:00:00
        #PBS -q checkpt
        #PBS -N pargenes_zip

        cd $PBS_O_WORKDIR

        tar -czf trimal-internal-odont-131-fasta-reduced-pargenes-bootstraps.tar.gz trimal-internal-odont-131-fasta-reduced-pargenes-bootstraps