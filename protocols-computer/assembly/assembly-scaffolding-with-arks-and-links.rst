.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _AssemblyScaffoldingWithArksAndLinks:

Assembly Scaffolding With Arks And Links
========================================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Assembly Scaffolding With Arks And Links`_.

.. _Assembly Scaffolding With Arks And Links: https://github.com/faircloth-lab/protocols/commits/master/protocols-computer/assembly/assembly-scaffolding-with-arks-and-links.rst


Purpose
-------

After assembling PacBio data with a program like canu_, we generally want to try and scaffold those contigs to achieve higher levels of assembly contiguity.  We can scaffold PacBio data using 10X linked reads.

Steps
-----

#. Compile arks_, links_, and install tigmint_ according to :ref:`CompilingArksAndLinks`. If you are in our lab, you don't need to do this - the directory is shared.

#. Before using the 10X linked read data that we have, we need to install the longranger_ program from 10X genomics.  That's a pretty easy process - you just need to go to the longranger_ website, grab the link and download that file to a reasonable location (in our case the shared directory ``/home/brant/project/shared/bin/``).

    .. code-block:: bash
    
        wget -O longranger-2.2.2.tar.gz "<long link from website>"
        tar -xzvf longranger-2.2.2.tar.gz

#. Setup a working directory.  Here, we're working with *Diglossa*, so:

    .. code-block:: bash

        mkdir 10x-diglossa-scaffold && cd $_


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

#. After running that, we need to generate a file of barcode multiplicities.  We can do that with a perl script from the arks_ package.  Before running this perl script, you need to create a configuration file containing the path to the processed linked read data from above.  In our ``10x-diglossa-scaffold`` directory, create a new directory for these data and create the ``reads.fof``:

    .. code-block:: text

      mkdir barcode-multiplicities && cd $_
      echo `readlink -f ../longranger/10x-diglossa/outs/barcoded.fastq.gz` > reads.fof

#. For example, my ``reads.fof`` will contain a single line that looks like:

    .. code-block:: text

        /ddnB/work/brant/10x-diglossa-scaffold/longranger-output/10x-diglossa/outs/barcoded.fastq.gz

#. And my overall directory structure will look like:

    .. code-block:: text

        .
        ├── barcode-multiplicities
        │   └── reads.fof
        └── longranger-output
            ├── 10x-diglossa
            ├── arks-make.orig
            ├── arks-make.txt
            ├── longranger_basic.err
            ├── longranger_basic.out
            ├── longranger.qsub
            └── raw-fastq

#. Now that we've created this file of filenames (FOFN or fofn), we can compute the barcode multiplicities.  This takes about 30 minutes for a 50 GB file of reads:

    .. code-block:: bash

        #!/bin/bash
        #PBS -q single
        #PBS -A <allocation>
        #PBS -l walltime=24:00:00
        #PBS -l nodes=1:ppn=1
        #PBS -V
        #PBS -N arks_multiplicities
        #PBS -o arks_multiplicities.out
        #PBS -e arks_multiplicities.err

        module load perl/5.24.0/INTEL-18.0.0

        export PATH=/home/brant/project/shared/bin/:$PATH

        cd $PBS_O_WORKDIR

        calcBarcodeMultiplicities.pl reads.fof > read_multiplicities.csv

#. Before we scaffold, we need to upload the contig files/pacbio assembly that we want to scaffold:

    .. code-block:: bash

        mkdir to-scaffold
        # rsync up to this directory from wherever contigs are located
        rsync -avLP diglossa.contigs.fasta brant@mike.hpc.lsu.edu:/home/brant/work/10x-diglossa-scaffold/to-scaffold/

#. So, now our directory structure looks something like:

    .. code-block:: bash

        .
        ├── barcode-multiplicities
        │   ├── multiplicities.qsub
        │   ├── read_multiplicities.csv
        │   └── reads.fof
        ├── longranger-output
        │   ├── 10x-diglossa
        │   ├── arks-make.orig
        │   ├── arks-make.txt
        │   ├── longranger_basic.err
        │   ├── longranger_basic.out
        │   ├── longranger.qsub
        │   └── raw-fastq
        └── to-scaffold
            └── diglossa.contigs.fasta

#. Finally, we are ready to run arks_.  Within your working directory, create a final a directory to hold the arks_ output:

    .. code-block:: bash

        mkdir arks-scaffolded && cd $_

#. arks_ is primarily run through a makefile, an example of which is `available on the arks github page <https://github.com/bcgsc/arks/blob/master/Examples/arks-make>`_.

#. I've already partially edited this make file to make arks run more easily given the way we have it installed.  You can download my edited version `here <https://gist.githubusercontent.com/brantfaircloth/a714928d2824a83684254587255f4c57/raw/de6050cf759b81e87342fc521a4cfb416401b333/arks-make.txt>`_.  We can see the parameters the makefile accepts by downloading the file and running it with ``make``:

    .. code-block:: bash

        wget -O arks-make.txt https://gist.githubusercontent.com/brantfaircloth/a714928d2824a83684254587255f4c57/raw/de6050cf759b81e87342fc521a4cfb416401b333/arks-make.txt
        make -f arks-make.txt

#. Generally speaking, there are some values in the makefile that we want to change, specifically the options for numbers of threads for both ``bwa`` and for ``arks``, which are named ``-t`` (around line 24) and ``-threads`` (around line 37). You want to adjust these values to the number of cores on whatever HPC system you are running on. arks_ doesn't use MPI, so you'll only submit to a single node (thus, you need to know how many processes that single node can run).

#. arks_ also does not do well with paths in its invocation, so in ``arks-scaffolded``, create symlinks to our fastq data and our assembly:

    .. code-block:: bash

        ln -s ../longranger-output/10x-diglossa/outs/barcoded.fastq.gz
        ln -s ../to-scaffold/diglossa.contigs.fasta

#. Now that we've edited the ``makefile``, we can setup the ``qsub`` file for the arks_ run, here assuming we're running on @supermike.  Of importance is the way that arks_ handles the expected file names - be sure to structure those correctly or you will get `this error <https://github.com/bcgsc/arks/issues/15>`_.  This essentially means that you need to refer to the contigs without including the ``.fasta`` extension and the reads without including the ``.fastq.gz`` extension:

    .. code-block:: bash

        #!/bin/bash

        #PBS -q checkpt
        #PBS -A <allocation>
        #PBS -l walltime=72:00:00
        #PBS -l nodes=1:ppn=16
        #PBS -V
        #PBS -N arks_scaffolding
        #PBS -o arks_scaffolding.out
        #PBS -e arks_scaffolding.err

        # load some modules
        module load gcc/6.4.0
        module load perl/5.24.0/INTEL-18.0.0
        module load boost/1.63.0/INTEL-18.0.0

        # activate the conda env in which we have installed tigmint
        conda activate scaffolding

        # make sure we inject all the correct paths
        export PATH=/home/brant/project/shared/bin/:/home/brant/project/shared/src/links_v1.8.7:$PATH

        cd $PBS_O_WORKDIR

        make -f arks-make.txt arks-tigmint \
            draft=diglossa.contigs \
            reads=barcoded \
            m=50-30000 o=3 time=1