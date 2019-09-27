.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running IQ-Tree
===============

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running IQ-Tree`_ 

.. _Running IQ-Tree: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-iq-tree.rst

Purpose
-------

IQ-Tree_ is a program to infer ML trees from multiple sequence alignment and to perform analyses of those (and other) trees.  We can analyze data in lots of different ways.  Here, we'll cover using IQ-Tree_ to analyze a concatenated data set and also to analyze a set of alignments to infer gene trees.


.. _IQ-Tree: http://www.iqtree.org

Preliminary Steps
-----------------

#. To compile Pargenes, see :ref:`CompilingIQtree`

Data Preparation
----------------

#. IQ-Tree accepts data in a multitude of formats: PHYLIP, FASTA, NEXUS, CLUSTALW.  Setup directory structure on @hpc to contain these data. Generally speaking, I make a project-specific folder within my ``work`` directory (where ``work`` is symlinked to ``/work/brant``.  So, using some fish data as an example:

    .. code-block:: bash

        mkdir work/fish-analyses
        cd fish-analyses
        mkdir alignments

#. On the transfer machine (@tabasco), navigate to the directory holding the alignment files and transfer the alignments files to @hpc (in this case, @supermic):

    .. code-block:: bash

        rsync -avLP ./ user@mike.hpc.lsu.edu:/home/brant/work/fish-analyses/alignments

.. admonition:: Note

        Unlike RAxML-NG, IQ-Tree does not require initial conversion of the data to a binary file.


Concatenated Alignments
-----------------------

This is pretty straightforward.  One thing to keep in mind is that IQ-Tree_ parallelizes poorly when you are analyzing partitioned alignments - these I would usually run in raxml-ng_.


Inferring the Best ML Tree (with ultra-fast bootstrapping)
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#. Sync up your concatenated alignment file.  IQ-Tree does not require that you turn your alignment into a binary.
#. Setup the following qsub file.  Note that we're simply appyling GTR+Gamma to the entire concatenated alignment here:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <your allocation>
        #PBS -l walltime=24:00:00
        #PBS -l nodes=13:ppn=20
        #PBS -q checkpt
        #PBS -V
        #PBS -N iqtree-nopart
        #PBS -o iqtree-nopart-mpi.out
        #PBS -e iqtree-nopart-mpi.err

        module load gcc/6.4.0
        module load impi/2018.0.128/intel64

        export TASKS_PER_HOST=1  # number of MPI tasks per host
        export THREADS_HOST=20   # number of threads spawned by each task on the host

        cd $PBS_O_WORKDIR

        mpirun -perhost ${NPERNODE:=1} -np ${PBS_NUM_NODES} -hostfile $PBS_NODEFILE \
            /project/brant/shared/bin/iqtree-omp-mpi \
            -s <your_alignment>.phy \
            -m GTR+G \
            -bb 1000 \
            -nt $THREADS_HOST



Individual Alignments
---------------------

For individual alignments, the data that we uploaded consist of a directory of aligmments (zipped or unzipped, doesn't really matter).  We basically need to run IQ-Tree_ against each of the alignment files in the directory.  As it does this, it will infer the best substitution model for each locus, then use that to infer the tree.

#. Prep a list of alignment files that we will feed to IQ-Tree_ so that it know which trees we want to infer.

    .. code-block:: bash

        # remove any existing list of loci to align
        rm trees-to-generate.list
        # loop over loci to build an input file
        for FULLPATH in $PWD/alignment-files-clean/*; do
            echo "$FULLPATH" >> trees-to-generate.list
        done

#. Prep the script that will actually run IQ-Tree_ named ``iqtree.sh``.  As written, this does not run any bootstrapping for a given locus.  If you want to do that, add the ``-bb 1000`` option to infer 1000 ultra-fast bootstraps.  Note that we are setting the number of cores needed for each alignment to 2 here (you may need to increase as needed for larger alignments):

    .. code-block:: bash

        #!/bin/bash

        ## set this manually
        CORES_PER_JOB=2

        ## DO NOT EDIT - this comes as input from GNU parallel on STDIN, via the sample.list file

        ALIGN=$1

        ## Here are the specific commands we are running
        # run iqtree
        /project/brant/shared/bin/iqtree -s $ALIGN -nt $CORES_PER_JOB

#. Make sure to make this script executable: ``chmod 0755 iqtree.sh``
#. Create a qsub script to run the job and specify the number of nodes/cores you will need.  Note that we are also declaring the number of CPUs per alignment here:

    .. code-block:: bash

        #PBS -A <your allocation>
        #PBS -l nodes=1:ppn=20
        #PBS -l walltime=2:00:00
        #PBS -q checkpt
        #PBS -N multi_iqtree

        module load gnuparallel/20170122

        # Number of Cores per job (needs to be multiple of 2)
        export CORES_PER_JOB=2

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --colsep '\,' \
                --progress \
                --joblog logfile.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a trees-to-generate.list \
                ./iqtree.sh {$1}

#. Submit the job: ``qsub pasta.qsub``.  Be sure to monitor the job with ``checkjob -j <job_number>`` to ensure you are using resources approriately.  You can also see how many trees have been inferred by running ``ls alignment-files-clean/*.treefile | wc -l``.


