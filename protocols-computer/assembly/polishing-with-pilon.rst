.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _PolishingWithPilon:

Polishing Assemblies with Pilon
===============================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Polishing Assemblies with Pilon`_.

.. _Polishing Assemblies with Pilon: https://github.com/faircloth-lab/protocols/commits/master/protocols-computer/assembly/polishing-with-pilon.rst


Purpose
-------

PacBio assemblies are great because they produce contigs with high contiguity.  However, the coverage of PacBio reads is something less than we desire, and PacBio reads are more prone error than something like Illumina reads.  We can use the Illumina reads to "polish" contigs produce from the "noisier" PacBio chemistry.  We can also use 10X_ genomics reads, if we have them, do to the same - because those are simply large insert Illumina reads.

.. note:: This protocol assumes you are using 10X_ reads.

.. warning:: The following assumes you are running on @qb2, which you should probably be doing because we need the ``bigmem`` queue for vertebrate-sized genomes.

Steps
-----

.. note:: You may have already performed some of these steps such as processing 10X_ reads to remove the barcode information or mapping the 10X_ reads to the assembly you'd like to polish.  Simply skip those steps if you have already performed them.

#. If you are using 10X_ reads, you need to process those to remove the internal barcodes and adapters.  That's most easily accomplished by installing and using longranger_. That's a pretty easy process - you just need to go to the longranger_ website, grab the link and download that file to a reasonable location (in our case the shared directory ``/home/brant/project/shared/bin/``, where it is already installed).

    .. code-block:: bash
    
        wget -O longranger-2.2.2.tar.gz "<long link from website>"
        tar -xzvf longranger-2.2.2.tar.gz

#. Setup a working directory:

    .. code-block:: bash

        mkdir pacbio-polish && cd $_


#. Now, we're going to run longranger_ to do some basic processing of the 10X linked read data.  Go ahead and make a directory within `10x-diglossa-scaffold` to hold the data:

    .. code-block:: bash

        mkdir longranger-ouput && cd $_


#. Prepare a ``qsub`` script to run longranger_ and process the reads.  This should take <24 hours for ~40 GB zipped sequence data.  The processing basically trims the reads to remove the barcode and adapter information and puts the barcode info in the fastq header:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q checkpt
        #PBS -A <allocation>
        #PBS -l walltime=24:00:00
        #PBS -l nodes=1:ppn=20
        #PBS -V
        #PBS -N longranger_basic
        #PBS -o longranger_basic.out
        #PBS -e longranger_basic.err
        #PBS -m abe
        #PBS -M brant@faircloth-lab.org

        export PATH=/home/brant/project/shared/bin/longranger-2.2.2/:$PATH

        cd $PBS_O_WORKDIR

        longranger basic \
            --id=<my_name> \
            --fastqs=/path/to/my/demuxed/fastq/files \
            --localcores 20 1>longranger-basic.stdout 2>longranger-basic.stderr

#. Once the reads have been processed, we want to map them to our genome assembly using ``bwa-mem`` and ``samtools``.  We can get all those installed (along with Pilon_) by creating a conda_ environment:

    .. code-block:: bash

        conda create -n polishing pilon bwa samtools

#. Now, we need to map the reads over.  The easiest thing to do is probably to create a new directory in ``pacbio-polished`` named ``bwa-aligned``

    .. code-block:: bash
        
        mkdir bwa-aligned && cd $_

#. Now, symlink in the `*.fastq.gz` file that we just created:

    .. code-block:: bash

        ln -s ../longranger-ouput/<my_name>/outs/barcoded.fastq.gz

#. Upload the assembly to this directory, using a tool like ``rsync``.  Here, we've uploaded ``diglossa.contigs.fa`` to the same directory (``bwa-aligned``) that contains our symlink to ``barcoded.fastq.gz``
#. Once that's uploaded, create a ``qsub`` script to run the ``bwa`` mapping job:


    .. code-block:: bash

        #!/bin/bash
        #PBS -q checkpt
        #PBS -A <allocation>
        #PBS -l walltime=36:00:00
        #PBS -l nodes=1:ppn=20
        #PBS -V
        #PBS -N bwa_mem
        #PBS -o bwa_mem.out
        #PBS -e bwa_mem.err

        source activate polishing

        cd $PBS_O_WORKDIR

        # index the assembly for bwa
        bwa index diglossa.contigs.fa

        # run bwa, use 20 threads for aligning and sorting and set memory for samtools at 3G per thread
        bwa mem -t 20 diglossa.contigs.fa barcoded.fq.gz | samtools sort -@20 -m 3G -o diglossa.contigs.barcoded.bam -
        samtools index diglossa.contigs.barcoded.bam


#. Submit that job and let it run.  It will take a fair amount of time (~18 hours to map something like 40 GB data)
#. Once the job has run, you should have a directory ``bwa-aligned`` that contains the output bam file, the reads, and the assembly.  Moving forward, we only need to care about the BAM file and the assembly.
#. Running Pilon_ is pretty simple - it just needs to use a lot of RAM (why we need to run it @qb2).  We need to setup an appropriate ``qsub`` script for the run:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q bigmem
        #PBS -A <allocation>
        #PBS -l walltime=72:00:00
        #PBS -l nodes=1:ppn=48
        #PBS -V
        #PBS -N pilon
        #PBS -o pilon.out
        #PBS -e pilon.err

        source activate polishing

        cd $PBS_O_WORKDIR

        # run pilon
        pilon -Xmx1400G --genome diglossa.contigs.fa \
            --bam diglossa.contigs.barcoded.bam \
            --changes --vcf --diploid --threads 48 \
            --output diglossa.contigs.polished
