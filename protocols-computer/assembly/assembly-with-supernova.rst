.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Assembly With Supernova
=======================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Assembly With Supernova`_ 

.. _Assembly With Supernova: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/assembly/assembly-with-supernova.rst

Purpose
-------

Supernova_ is a program for assembling 10X Genomics Linked Read Data.


.. _Supernova: https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/installation

Preliminary Steps
-----------------

#. To install Supernova, see :ref:`InstallingSupernova`

Steps
-----

#. Prior to running Supernova_, it's a good idea to get an idea of the count of reads that you have for a given sample.  You want to be inputting roughly 56-60X coverage, per the 10X instructions.  You can compute the counts of reads that you have using:

    .. code-block:: bash

        for i in clean-reads/*; do echo $i; gunzip -c $i/split-adapter-quality-trimmed/*-READ1.fastq.gz | wc -l | awk '{print $1/4}'; done

#. This will output a count of R1 reads by sample to the console.  To get the total counts of reads, multiple by 2.  To get a rough estimate of coverage, multiply that by the length of both reads.  Divide that number by the size of your genome to get some idea of coverage.  We can dial down the number of reads when we run Supernova_ if we need to.  Guidance regarding the number of reads to use with Supernova_ can be found at `this page <https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/running>`_.

#. Setup a submission script for QB2 (in our case).  Generally speaking, avian-sized genome assemblies are going to need something like 256 GB of RAM, whereas mammal sized genomes may need up to 512.  However, Supernova should be run on **AT LEAST** 16 CPU cores, and we want it to finish in a reasonable amount of time (< 72 hours).  So, on QB2, that means we'll run a job with 18 of the 48 cores available on a QB2 bigmem node.  This will net us ~562 GB RAM.  Because of the way the program runs, we need to explicitly limit the number of cores and RAM used by the Supernova_ process. We'll slightly undershoot the total RAM allocated to the job (limiting it to 512 GB of the 562 GB).

    .. code-block:: bash

        #!/bin/bash
        #PBS -q bigmem
        #PBS -A <allocation>
        #PBS -l walltime=02:00:00
        #PBS -l nodes=1:ppn=18
        #PBS -V
        #PBS -N supernova_assembly
        #PBS -o supernova_assembly.out
        #PBS -e supernova_assembly.err


        export PATH=$HOME/bin/supernova-2.1.1:$PATH

        cd $PBS_O_WORKDIR
        supernova run \
            --id=<my_assembly_name> \
            --fastqs=/path/to/my/demuxed/fastq/files \
            --maxreads=<maxreads determined based on above> \
            --localcores 18 \
            --localmem 512
