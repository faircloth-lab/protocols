.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running RAxML-NG
================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running RAxML-NG`_ 

.. _Running RAxML-NG: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-raxml-ng.rst

Purpose
-------

We use RAxML-NG (formerly RAxML or ExaML) to infer maximum likelihood trees from multiple sequence alignment data.  RAxML-NG is now one of the preferred options because it can run across multiple nodes and cores and also take, as input, a number of evolutionary models estimated by a program like `modeltest-ng`_

.. _modeltest-ng: https://github.com/ddarriba/modeltest


Citation
--------

`RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference <http://doi.org/10.1093/bioinformatics/btz305>`_


Preliminary Steps
-----------------

#. To compile RAxML-NG, see :ref:`CompilingRaxmlng`


Data Preparation
----------------

#. RAxML accepts data in two format: PHYLIP and FASTA. Setup directory structure on @supermike to contain these data. Generally speaking, I make a ``project`` folder within my ``work`` directory (where ``work`` is symlinked to ``/work/brant``.  So, using Anna's Diglossa as an example:

    .. code-block:: bash

        mkdir work/anna-diglossa
        cd anna-diglossa
        mkdir alignments

#. On the transfer machine (@tabasco), navgate to the directory holding the alignment files and transfer the alignments files to @supermike:

    .. code-block:: bash

        rsync -avLP ./ user@mike.hpc.lsu.edu:/home/brant/work/anna-diglossa/alignments

#. Now that that's finished, setup a file that will (1) convert the alignment to a binary format, and (2) estimate the number of nodes/cores needed for optimal analysis. I do this in a file named something like ``raxml-parse.qsub``:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q single
        #PBS -A <allocation>
        #PBS -l walltime=02:00:00
        #PBS -l nodes=1:ppn=1
        #PBS -V
        #PBS -N raxmlng-parse
        #PBS -o raxmlng-parse.out
        #PBS -e raxmlng-parse.err

        module load gcc/6.4.0
        module load impi/2018.0.128

        cd $PBS_O_WORKDIR

        /project/brant/shared/bin/raxml-ng \
        --msa /path/to/alignment/alignment.phylip \
        --model GTR+G \
        --parse

    .. admonition:: Note

        RAxML-NG is different from earlier versions of RAxML because you now have the ability to specify many, many different models of sequence evolution.  The options for evolution models are `on the RAxML-NG website <https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model>`_ and should be perused. It's also possible to specify multiple models using a partition file (``partition.txt``) that looks something like:

        .. code-block:: text

            JC+G, p1 = 1-100, 252-400
            HKY+F, p2 = 101-180, 251
            GTR+I, p3 = 181-250

        And then creating the binary alignment file with a command similar to:

        .. code-block:: bash

            /project/brant/shared/bin/raxml-ng \
            --msa /path/to/alignment/alignment.phylip \
            --model partition.txt \
            --parse

        When creating the binary alignment file, you specify the model to use for the given data set.  This model will be carried over to all subsequent analyses using this binary alignment file - which is why we don't specify particular models in the sections below.

#. This will produce binary alignment files within ``/path/to/alignment/``.  These files will have an ``.rba`` extension (so the file created here was ``drop2-mafft-nexus-edge-trimmed-clean-75p.phylip.raxml.rba``).  To look at other information regarding the aligment (particularly how many nodes/cores to use), open up ``raxmlng-parse.out`` with something like ``less``.  You should see content that looks like:

    .. code-block:: text

        Analysis options:
        run mode: Alignment parsing and compression
        start tree(s):
        random seed: 1558381540
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX
        parallelization: PTHREADS (8 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: alignments/drop2-mafft-nexus-edge-trimmed-clean-75p.phylip
        [00:00:00] Loaded alignment with 116 taxa and 2094052 sites

        WARNING: Fully undetermined columns found: 31736

        NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed)
        NOTE: was saved to: /ddnB/work/brant/anna-diglossa/alignments/drop2-mafft-nexus-edge-trimmed-clean-75p.phylip.raxml.reduced.phy

        Alignment comprises 1 partitions and 720183 patterns

        Partition 0: noname
        Model: GTR+FO+G4m
        Alignment sites / patterns: 2062316 / 720183
        Gaps: 25.14 %
        Invariant sites: 88.39 %


        NOTE: Binary MSA file created: alignments/drop2-mafft-nexus-edge-trimmed-clean-75p.phylip.raxml.rba

        * Estimated memory requirements                : 20220 MB

        * Recommended number of threads / MPI processes: 96

        Please note that numbers given above are rough estimates only.
        Actual memory consumption and parallel performance on your system may differ!

Inferring the Best ML Tree (with bootstrapping)
-----------------------------------------------

After creating the binary alignment file and getting an idea of the number of MPI processes that are needed, you need to infer the tree, ideally with some support values.

You have several ways of doing this, one of which is to use what I call "standard" MPI mode, which just gives RAxML a number of CPUs to spread the data across, and all the CPUs talk to each other over the interconnects using MPI.

.. warning::

    **Still testing.**

    The other way of setting up the run is to use what's known as "hybrid" mode, which combines parallelization across HPC nodes (using MPI) with parallelization within nodes (using Pthreads).


.. admonition:: Note

    In all of the following, you can prefix the name of your output files by adding the argument ``--prefix <some name>``.  And, when generating consensus trees, you can root those on some outgroup using ``--outgroup taxon1,taxon2,taxon3,..., taxonQQQ``.


Standard MPI Mode
:::::::::::::::::

Using Standard MPI Mode To Search for the Best ML Tree + Boostrapping
.....................................................................

Given the core count and RAM usage estimated above, on @supermike, we need to run 6 nodes each with 16 CPUS for a total of 96 CPUs.  We will also set this run up to automatically search for both the *best* ML tree and **bootstrap replicates** for this best ML tree.  That's accomplished with the ``--all`` option.  The other option we are passing is the ``--best-trees autoMRE`` option, which will generate bootstrap trees until those converge.  If you need to set the maximum number of boostrap replicated to generate using autoMRE, specify ``--bs-trees autoMRE{500}``, which will limit the analyses to only 500 trees (default is 1000). The ``--seed`` that we're setting (which we pass as an environment variable ``$SEED`` whose value it taken from $RANDOM) let's us repeat the exact analysis in the future, if needed.

With that information in hand, setup a second submission script ``raxml-best-tree.qsub`` that contains a version of the following:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q checkpt
        #PBS -A <your_allocation>
        #PBS -l walltime=72:00:00
        #PBS -l nodes=6:ppn=16
        #PBS -V
        #PBS -N raxmlng-std-mpi
        #PBS -o raxmlng-std-mpi.out
        #PBS -e raxmlng-std-mpi.err

        module load gcc/6.4.0
        module load impi/2018.0.128

        cd $PBS_O_WORKDIR
        SEED=$RANDOM
        echo $SEED

        mpiexec -np 96 -machinefile $PBS_NODEFILE /project/brant/shared/bin/raxml-ng-mpi \
            --msa alignments/drop2-mafft-nexus-edge-trimmed-clean-75p.phylip.raxml.rba \
            --seed $SEED \
            --all \
            --bs-trees autoMRE

Using Standard MPI Mode To Search for the Best ML Tree
......................................................

Sometimes, the tree you are trying to infer is large (due to the # of tips, the amount of data, or both), and you want to separate the inference of the best ML tree from the generation of bootstrap replicates.  To infer only the best ML tree, use something like the following.  The ``--search`` option tells RAxML to do the tree search and ``--seed`` is described above.

    .. code-block:: bash

        #!/bin/bash
        #PBS -q checkpt
        #PBS -A <your_allocation>
        #PBS -l walltime=72:00:00
        #PBS -l nodes=6:ppn=16
        #PBS -V
        #PBS -N raxmlng-std-mpi
        #PBS -o raxmlng-std-mpi.out
        #PBS -e raxmlng-std-mpi.err

        module load gcc/6.4.0
        module load impi/2018.0.128

        cd $PBS_O_WORKDIR
        SEED=$RANDOM
        echo $SEED

        mpiexec -np 96 -machinefile $PBS_NODEFILE /project/brant/shared/bin/raxml-ng-mpi \
            --msa alignments/drop2-mafft-nexus-edge-trimmed-clean-75p.phylip.raxml.rba \
            --seed $SEED \
            --search



Using Standard MPI Mode To Bootstrap
....................................

Along similar lines, if you've separated how RAxML runs into two parts, you would run the boostrapping for a particular set of data using the following.  The ``--seed`` argument is described above, the ``--bootstrap`` argument tells RAxML to do bootstrapping, the the ``--bs-trees`` argument is described above:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q checkpt
        #PBS -A <your_allocation>
        #PBS -l walltime=72:00:00
        #PBS -l nodes=6:ppn=16
        #PBS -V
        #PBS -N raxmlng-std-mpi
        #PBS -o raxmlng-std-mpi.out
        #PBS -e raxmlng-std-mpi.err

        module load gcc/6.4.0
        module load impi/2018.0.128

        cd $PBS_O_WORKDIR
        SEED=$RANDOM
        echo $SEED

        mpiexec -np 96 -machinefile $PBS_NODEFILE /project/brant/shared/bin/raxml-ng-mpi \
            --msa alignments/drop2-mafft-nexus-edge-trimmed-clean-75p.phylip.raxml.rba \
            --seed $SEED \
            --bootstrap \
            --bs-trees autoMRE

Integrating the Best ML Tree with the Bootstraps
................................................

And, if you have separate files for the best ML tree and the boostrap replicates, you can integrate those using (**Note** that I'm using the ``single`` queue here with a very short walltime because this runs quickly).


    .. code-block:: bash

        #!/bin/bash
        #PBS -q single
        #PBS -A <your_allocation>
        #PBS -l walltime=00:10:00
        #PBS -l nodes=1:ppn=1
        #PBS -V
        #PBS -N raxmlng-std-mpi
        #PBS -o raxmlng-std-mpi.out
        #PBS -e raxmlng-std-mpi.err

        module load gcc/6.4.0
        module load impi/2018.0.128

        cd $PBS_O_WORKDIR

        /project/brant/shared/bin/raxml-ng \
        --tree /path/to/bestML.tree \
        --bs-trees /path/to/bootstraps.tree \
        --support


Post-hoc Tree Evaluation
------------------------

Sometimes after the best ML tree search, you will see ML trees inferred with different likelihood values.  It may be important to evaluate the differences among the best ML trees.  It may also be important for you to do things like compute concordance factors when comparing results from concatenated trees to something like gene tree topologies.

Examining Likelihoods and RF Distance
:::::::::::::::::::::::::::::::::::::

You can easily compute the RF distance among the best ML trees inferred using RAxML

