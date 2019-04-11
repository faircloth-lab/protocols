.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Assembly With Itero
===================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Assembly With Itero`_ 

.. _Assembly With Itero: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/assembly/assembly-with-itero.rst

Purpose
-------

itero_ is a pipeline to generate gene trees from a large set of loci, using the most appropriate site rate substitution model.


.. _itero: https://github.com/faircloth-lab/itero
.. _manual: http://itero.readthedocs.io

Preliminary Steps
-----------------

#. To install itero_, please see the manual_.  If you plan to use a cluster/HPC, then be sure to install itero_ there.

Steps
-----

#. Prior to running itero_, it's a good idea to get an idea of the count of reads that you have for a given sample.  If this number is > 3M reads per sample and these are UCE data, we probably want to downsample the data to a manageable size - something like 3M total reads (1.5 M R1 and 1.5M R2 reads) per sample.  You also probably want to do this **after trimming**. Once you've cleaned the reads, you can compute the count of reads by:

    .. code-block:: bash

        for i in clean-reads/*; do echo $i; gunzip -c $i/split-adapter-quality-trimmed/*-READ1.fastq.gz | wc -l | awk '{print $1/4}'; done

#. This will output a count of R1 reads by sample to the console, and you can use a regular expression to re-arrange those bits into a CSV file.  Load that CSV file in excel (or sort in some manner) so that you can determine which samples have 3M reads (1.5M read each for ``R1`` and ``R2``).

#. One you have the list of those samples you'd like to downsample, you need to create a text file to hold their names, call it something like ``samples-to-downsample.txt``:

    .. code-block:: text

        alectura-lathami2
        anas-platyrhynchos
        anser-erythropus
        anseranas-semipalmata
        biziura-lobata
        chauna-torquata
        colinus-cristatus
        coturnix-coturnix
        crax-alector
        gallus-gallus
        malacorhynchus-membranaceus
        megapodius-eremita
        numida-meleagris
        oxyura-jamaicensis
        rollulus-rouloul

#. In your working directory, imagine you have ``clean-reads`` containing your trimmed read data.  Create a new directory ``downsampled-reads``, and ``cd`` into that.  Make sure the text file from above is in this new directory. Now, assuming you have ``seqtk`` in your ``$PATH`` somewhere:

    .. code-block:: bash

        for sample in ``cat samples-to-downsample.txt``;
            do rnum=$RANDOM;
            reads=1500000;
            echo 'sampling: ' ${sample} ${reads};
            echo 'using: ' ${rnum};
            mkdir ${sample};
            seqtk sample -s $rnum ../clean-reads/${sample}/split-adapter-quality-trimmed/${sample}-READ1.fastq.gz $reads | gzip > ./${sample}/${sample}-READ1.${reads}.fastq.gz;
            seqtk sample -s $rnum ../clean-reads/${sample}/split-adapter-quality-trimmed/${sample}-READ2.fastq.gz $reads | gzip > ./${sample}/${sample}-READ2.${reads}.fastq.gz;
        done

#. This will create files containing 1,500,000 reads that have been randomly sampled from your larger population of reads.  If you need to, symlink in ``R1`` and ``R2`` files from samples having fewer than 3M total reads.  And, if you need to, upload all of those reads to wherever you are running your assembly (e.g. ``supermic``).

#. On ``supermic``, you generally want to split your samples up into batches of about 20 taxa, and you want to run ~2 batches of 20 taxa at the same time (itero_ uses a lot of IO, so we want to be nice and not suck all the available IO up).  Based on your read locations, you want to create a configuration file for a batch of assemblies.  That looks something like this:

    .. code-block:: text

        [reference]
        /home/brant/work/eb2/uce-5k-probes.loci.fasta

        [individuals]
        alectura-lathami2:/home/brant/work/eb2/batch-1/raw-reads/alectura-lathami2
        anas-platyrhynchos:/home/brant/work/eb2/batch-1/raw-reads/anas-platyrhynchos
        anser-erythropus:/home/brant/work/eb2/batch-1/raw-reads/anser-erythropus
        anseranas-semipalmata:/home/brant/work/eb2/batch-1/raw-reads/anseranas-semipalmata
        biziura-lobata:/home/brant/work/eb2/batch-1/raw-reads/biziura-lobata
        chauna-torquata:/home/brant/work/eb2/batch-1/raw-reads/chauna-torquata
        colinus-cristatus:/home/brant/work/eb2/batch-1/raw-reads/colinus-cristatus
        coturnix-coturnix:/home/brant/work/eb2/batch-1/raw-reads/coturnix-coturnix
        crax-alector:/home/brant/work/eb2/batch-1/raw-reads/crax-alector
        gallus-gallus:/home/brant/work/eb2/batch-1/raw-reads/gallus-gallus
        malacorhynchus-membranaceus:/home/brant/work/eb2/batch-1/raw-reads/malacorhynchus-membranaceus
        megapodius-eremita:/home/brant/work/eb2/batch-1/raw-reads/megapodius-eremita
        numida-meleagris:/home/brant/work/eb2/batch-1/raw-reads/numida-meleagris
        oxyura-jamaicensis:/home/brant/work/eb2/batch-1/raw-reads/oxyura-jamaicensis
        rollulus-rouloul:/home/brant/work/eb2/batch-1/raw-reads/rollulus-rouloul


#. Then, you need so create a job submission script that looks something like:

    .. code-block:: bash

        #PBS -A <allocation_name>
        #PBS -l nodes=5:ppn=20
        #PBS -l walltime=72:00:00
        #PBS -q checkpt
        #PBS -N itero_batch1

        ulimit -n 10000

        # move into the directory containing this script
        cd $PBS_O_WORKDIR
        echo $PBS_NODEFILE

        source activate itero
        mpirun -hostfile $PBS_NODEFILE -n 100 itero assemble mpi --config batch-1.conf --output batch-1-assembly --local-cores 20 --clean
        source deactivate itero

#. This will run the important bits using 100 CPUs (20 CPUs for the bwa steps).

#. If you saved the above submission script as ``itero.pbs``, then submit the job:

    .. code-block:: text

        qsub itero.pbs

#. You can also submit additional batches (beyond 2) and make them dependent on the earlier batches finishing - in this way you can submit lots of (batch) jobs, but only have 2 running at the same time and not hog resources.  If your submission script is called ``itero.pbs``, then you need to determine the ``job_id`` you want the new job to start after and run:

    .. code-block:: text

        qsub -W depend=afterok:<job_id> itero.pbs