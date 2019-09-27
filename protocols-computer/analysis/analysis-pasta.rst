.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running Pasta
=============

:Author: Carl Oliveros and Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running Pasta`_ 

.. _Running Pasta: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-pasta.rst

Purpose
-------

Pasta_ is essentially an updated version of Saté_, and Pasta_ may or may not be incorporated with Treeshrink_ (we strictly run Pasta_ here).  Pasta_ is an iterative aligner that tends to align better than simply using mafft_ on its own.

.. _Pasta: https://github.com/smirarab/pasta
.. _Saté: http://phylo.bio.ku.edu/software/sate/sate.html
.. _Treeshrink: https://github.com/uym2/TreeShrink

Preliminary Steps
-----------------

#. Installing Pasta_ is a little tricky.  If you have not already, create a ``conda`` environment for Pasta_

    .. code-block:: bash

        conda create -n pasta python=3

#. Switch to that environment, make a ``tmp`` directory in it, and pull down the pasta code and binaries:

    .. code-block:: bash

        source activate pasta
        cd ~/anaconda/envs/pasta/
        mkdir -p tmp/pasta-code
        cd tmp/pasta-code
        git clone https://github.com/smirarab/pasta.git
        git clone https://github.com/smirarab/sate-tools-linux.git
        cd pasta
        python setup.py develop

#. To leave the environment, run:

    .. code-block:: bash

        source deactivate pasta

Steps
-----

#. Pasta_ enters the equation when we're trying to align DNA sequences.  Here, I'll discuss these steps in the context of using the phyluce_ pipeline, although many of the steps in the approach are similar regardless of whether you are using phyluce_ or not.

#. You can implement Pasta_ at several stages of the phyluce_ pipeline, but the easiest is probably after you have identified the UCE loci and extracted those loci to what we call a "monolithic" fasta file.  First thing you want to do it "explode" that monolithic FASTA by locus:

    .. code-block:: bash

        phyluce_assembly_explode_get_fastas_file \
            --input my-monolithic.fasta \
            --output exploded-loci

#. Now what you should have is a directory of fasta files, one for each locus in your data set. You likely want to filter these loci to remove really short stuff - typically something like those sequences having < 50% of the median length of all sequences for a particular locus:

    .. code-block:: bash

        phyluce_assembly_filter_seqs_from_fastas \
            --input exploded-loci \
            --output exploded-loci-length-filtered \
            --filtered-sequences-file exploded-loci_fasta.shorts \
            --proportion 0.5 \
            --cores 12 \
            --log-path log

#. And, we want to filter those FASTA files to remove loci that have fewer than 4 taxa.  We can do this using some code meant for alignment data:

    .. code-block:: bash

        phyluce_align_filter_alignments \
            --alignments exploded-loci-length-filtered \
            --output exploded-loci-length-min-4-taxa-filtered \
            --min-taxa 4

#. Once that's done, we want to package up those alignments and sync them to one of the HPC clusters.

    .. code-block:: bash

        # package them up
        tar -czvf exploded-loci-length-min-4-taxa-filtered.tar.gz exploded-loci-length-min-4-taxa-filtered

        # sync to supermic
        rsync -avLP exploded-loci-length-min-4-taxa-filtered.tar.gz  you@smic.hpc.lsu.edu:/home/you/work/

#. On the HPC machine, unarchive those files:

    .. code-block:: bash

        tar -xzvf exploded-loci-length-min-4-taxa-filtered.tar.gz

#. Let's rename that directory, to make things simpler

    .. code-block:: bash

        mv exploded-loci-length-min-4-taxa-filtered exploded-loci

#. Because of the way that we need to parallelize Pasta_, we need to create a CSV file that contains the list of files we want to align, along with the output directory information, and the locus name.  You can most easily do this with a script like the following:

    .. code-block:: bash

        # choose a name for the output directory
        OUTDIR="pasta-output"
        # remove any existing list of loci to align
        rm loci-to-align.list
        # loop over loci to build an input file
        for FULLPATH in $PWD/exploded-loci/*; do
            FILE=`basename $FULLPATH`;
            NAME="${FILE%%.*}"
            echo "$FULLPATH,$OUTDIR,$NAME" >> loci-to-align.list
        done

#. This creates a file, ``loci-to-align.list`` that looks like the following:

    .. code-block:: text

        /home/brant/work/jarvis-align/exploded-loci/uce-1003.unaligned.fasta,pasta-output,uce-1003
        /home/brant/work/jarvis-align/exploded-loci/uce-1004.unaligned.fasta,pasta-output,uce-1004
        /home/brant/work/jarvis-align/exploded-loci/uce-1005.unaligned.fasta,pasta-output,uce-1005
        ...

#. Now, create a bash script named ``pasta.sh`` that GNU Parallel will call to run the individual alignments.  Note that we are setting the number of cores needed for each alignment to 2 here.  You may need to adjust this value if you have very many taxa in each alignment or very many alignments (to increase the amount of RAM per alignment).  We are also using the ``ginsi`` aligner from mafft_, which seems to deal with abberrant sections of sequence pretty well:

    .. code-block:: bash

        #!/bin/bash

        source activate pasta

        ## set this manually
        CORES_PER_JOB=2

        ## DO NOT EDIT - this comes as input from GNU parallel on STDIN, via the sample.list file

        FASTA=$1
        OUTDIR=$2/$3
        mkdir -p $OUTDIR

        ## Here are the specific commands we are running
        # run pasta
        python ~/anaconda/envs/pasta/bin/pasta/run_pasta.py --datatype=DNA --num-cpus=$CORES_PER_JOB --input=$FASTA --output-directory=$OUTDIR --aligner=ginsi

#. Make sure to make this script executable: ``chmod 0755 pasta.sh``
#. Create a qsub script to run the job and specify the number of nodes/cores you will need.  Note that we are also declaring the number of CPUs per alignment here:

    .. code-block:: bash

        #PBS -A <allocation_name>
        #PBS -l nodes=10:ppn=20
        #PBS -l walltime=12:00:00
        #PBS -q checkpt
        #PBS -N multi_pasta

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
                -a loci-to-align.list \
                ./pasta.sh {$1} {$2} {$3}

#. Submit the job: ``qsub pasta.qsub``.  Be sure to monitor the job with ``checkjob -j <job_number>`` to ensure you are using resources apprpriately.

#. One the jobs is finished, your alignment data should be in the ``pasta-output`` folder.  However, the way that Pasta_ formats things, they data are sort of messy.  So, you have a couple of options.  The easiest is probably to create a new folder and symlink the alignments for each loci to the new folder.  This is easiest to do if you are using zsh_:

    .. code-block:: bash

        mkdir alignments && cd alignments

        # BASH version (use if you are doing this on @hpc)
        for i in ../pasta-output/*; do 
            LOCUS=`basename $i`
            ln -s $i/*.aln $LOCUS.align.fasta;
        done

        # ZSH version (use if you are doing this on our analysis machines and you use ZSH)
        for i in ../pasta-output/*; do 
            ln -s $i/*.aln $i:t.align.fasta;
        done

#. Once we have our symlinks, we can make a copy of the ``alignments`` folder and follow the symlinks using a special command.  This will basically create a new directory containing correctly renamed alignments.

    .. code-block:: bash

        # must use the -rL option to copy links as real files
        cp -rL alignments alignment-files

#. Copy that back to whatever of our local machines you are using.

    .. code-block:: bash

        tar -czvf aligment-files.tar.gz alignment-files
        rsync -avLP you@smic.hpc.lsu.edu:/home/you/aligment-files.tar.gz ./

#. Finally, you will most likely want to trim the resulting alignments to remove aberrant sequences that ``ginsi`` has offset from the "good" parts of the alignment.

    .. code-block:: bash

        phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
            --alignments alignment-files \
            --output alignment-files-trim \
            --input-format fasta \
            --output-format nexus \
            --cores 12

#. Don't forget to strip the locus name information from each alignment:

    .. code-block:: bash

        phyluce_align_remove_locus_name_from_nexus_lines \
            --alignments alignment-files-trim \
            --output alignment-files-clean \
            --input-format nexus \
            --output-format nexus \
            --cores 12