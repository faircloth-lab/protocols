.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running RADcap Analysis
=======================

:Author: Jessie Salter and Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running RADcap Analysis`_ 

.. _Running RADcap Analysis: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-radcap.rst

Purpose
-------

The following assumes you are demultiplexing RADcap data prepared with enzymes and the i5-8N tag.  Otherwise, if you used standard libraries for RADcap locus enrichment, you can demultiplex those data like usual.

Preliminary Steps
-----------------

#. To compile stacks when using LSU HPC, be sure to `module load gcc/6.4.0` or enable this modules in your `~/.modules` file 
#. Get stacks, configure (w/ home directory install), and install.  The commands below **need to be modified** because they are set to install everything into ``/project/brant/home/``, which you don't have access to.  E.g.,

    .. code-block:: bash

        wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.54.tar.gz
        tar -cxzvf stacks-2.54.tar.gz
        cd stacks-2.54

        export CC=`which gcc`
        export CXX=`which g++`
        ./configure --prefix=/project/brant/home/
        make
        # if using an entire node, you can `make -j 20`
        make install

#. Get BBmap, and install that somewhere.  Basically download and place the files somewhere in your ``$PATH``

    .. code-block:: bash

        mkdir $HOME/bin
        wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.87.tar.gz
        tar -xzvf BBMap_38.87.tar.gz
        # this will create a folder bbmap which you need to add to your $PATH

.. admonition:: Note

        If you are in my lab group, these are installed in ``$HOME/project/brant/bin``


Steps
-----

#. Upload the relevant data to some location on @smic.  These should not have been demultiplexed in any way.

#. You may want to check to ensure the MD5 checksums of your files uploaded match the MD5 checksums that you expect.  Usually, you receive these from the sequencing center.

#. If you have multiple files (for some reason), you can combine the files together for READ1 and then combine the files together for READ2.

Your Data Contain Randomly Sheared DNA ("standard" libraries)
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#. If your data contain RADcap performed on randomly sheared or "standard" sequencing libraries that are mixed with "regular" RAD-cap libraries, we'll go ahead and demultiplex the randomly sheared, RADcap data, first. Once that's done, we will demultiplex the remaining reads containing i5-8N tags.

#. I am starting with a directory structure that looks like this:

    .. code-block:: bash

        .
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_1.fq.gz
        └── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_2.fq.gz


#. The procedure for these standard libraries is identical to the one described in :ref:`Demultiplexing a Sequencing Run`, so refer to that document, demultiplex, rename your files, and return here.

#. When I'm done with this first step, my directory structure looks something like this:

    .. code-block:: bash

        .
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_1.fq.gz
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_2.fq.gz
        └── random-libraries
            ├── ACCATCCA+ACATTGCG_R1_001.fastq.gz
            ├── ACCATCCA+ACATTGCG_R2_001.fastq.gz
            ├── ...
            ├── demuxbyname.e631878
            ├── demuxbyname.o631878
            ├── demux.qsub
            ├── my_barcodes.txt
            ├── Undetermined_R1_001.fastq.gz
            └── Undetermined_R2_001.fastq.gz


Your Data Contain i7 + i5-8N Tags
:::::::::::::::::::::::::::::::::

#. We next need to demultiplex the data that are ``i7 + i5-8N`` by the i7 tag, then we will process those "plates" worth of data, separately to find restriction sites, deal with the i5-8N tags, etc.  We can do this several ways and we can even use stacks, however stacks is relatively slow for this task when faster options exist, so we will (as above) use `demuxbyname.sh` from BBMap_.

#. Before we demultiplex, we need to create a file of indexes.  This will be similar to, but slightly different from how these data were demultiplexed, above.  So, create a file, e.g., `my_i7_indexes.txt` that contains all of the i7 indexes you have paired with i5-8N indexes.  The file **MUST** have each entry appended with a `+` because we are using substring matching to ensure we get what we want, and this makes the most appropriate substring.  So, if I have 4 i7 indexes, I want to create a directory containing a file that looks like the following, where the sequences to the left of the `+` represent the i7 indexes.  So, ``mkdir i5-8N_libraries``, then create a file `my_i7_indexes.txt` in that directory containing:

    .. code-block:: text

        CGAACTGT+
        CATTCGGT+
        TCGGTTAC+
        AGTCGCTT+

#. My directory structure now looks like:

    .. code-block:: text

        .
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_1.fq.gz
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_2.fq.gz
        ├── i5-8N_libraries
        │   └── my_i7_indexes.txt
        └── random-libraries
            ├── ACCATCCA+ACATTGCG_R1_001.fastq.gz
            ├── ACCATCCA+ACATTGCG_R2_001.fastq.gz
            ├── ...
            ├── demuxbyname.e631878
            ├── demuxbyname.o631878
            ├── demux.qsub
            ├── my_barcodes.txt
            ├── Undetermined_R1_001.fastq.gz
            └── Undetermined_R2_001.fastq.gz

#. We're almost ready to demultiplex, but before we do and if you already demultiplexed "standard" libraries, **make sure the files you are demultiplexing this time are the ``Undetermined_*`` files left over from the initial round of demultiplexing** (e.g. `Undetermined_R1_001.fastq.gz` in the directory structure, above).  Once you are sure that is so, you can setup demultiplexing.  I usually do this in a folder one level below the read data I am demultiplexing (``i5-8N_libraries``):

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <allocation>
        #PBS -l nodes=1:ppn=20
        #PBS -l walltime=12:00:00
        #PBS -q workq
        #PBS -N demuxbyname

        module load jdk/1.8.0_161

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        $HOME/project/shared/bin/bbmap/demuxbyname.sh \
            prefixmode=f \
            substring=t \
            in=../random-libraries/Undetermined_R1_001.fastq.gz \
            in2=../random-libraries/Undetermined_R2_001.fastq.gz \
            out=%_R1_001.fastq.gz \
            out2=%_R2_001.fastq.gz \
            outu=Undetermined_R1_001.fastq.gz \
            outu2=Undetermined_R2_001.fastq.gz \
            names=my_i7_indexes.txt

#. Once that finishes, the directory structure looks something like this (``random-libraries`` is collapsed):

    .. code-block:: text

        .
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_1.fq.gz
        ├── AEM1_CKDL200166465-1a_HF5GHCCX2_L3_2.fq.gz
        ├── i5-8N_libraries
        │   ├── AGTCGCTT+_R1_001.fastq.gz
        │   ├── AGTCGCTT+_R2_001.fastq.gz
        │   ├── CATTCGGT+_R1_001.fastq.gz
        │   ├── CATTCGGT+_R2_001.fastq.gz
        │   ├── CGAACTGT+_R1_001.fastq.gz
        │   ├── CGAACTGT+_R2_001.fastq.gz
        │   ├── demuxbyname.e631893
        │   ├── demuxbyname.o631893
        │   ├── demux.qsub
        │   ├── my_i7_indexes.txt
        │   ├── TCGGTTAC+_R1_001.fastq.gz
        │   ├── TCGGTTAC+_R2_001.fastq.gz
        │   ├── Undetermined_R1_001.fastq.gz
        │   └── Undetermined_R2_001.fastq.gz
        └── random-libraries

#. Now, within the ``i5-8N_libraries`` directory, we need to demultiplex the samples within each "plate" (or for each set of unique i7 indexes).  There are several ways that you can set this up (serially or parallel), but let's assume that you just want to do this serially for each "plate", because there are only a handful of plates.  I would start by making a directory for each plate (e.g. ``mkdir CGAACTGT+`` or rename as needed).  Once that's done, you need to navigate into that directory and create a file of the **internal** index sequences and the sample names associated with those sequences (e.g.``CGAACTGT_internal_indexes.txt``) .  You also need to add some nucleotides to each index, and this can get a little messy, to the easiest way to do this is to download `this file <https://www.dropbox.com/s/lr6p7k5894wghux/stacks-worksheet.xlsx?dl=1>`_, fill in the required information, and copy the contents noted into the ``CGAACTGT_internal_indexes.txt`` file.  The file contents will look something like:

    .. code-block:: text

        CCGAATG	CTAACGT	Sample_22
        CCGAATG	TCGGTACT	Sample_23
        CCGAATG	GATCGTTGT	Sample_24
        CCGAATG	AGCTACACTT	Sample_25
        CCGAATG	ACGCATT	Sample_26
        CCGAATG	GTATGCAT	Sample_27
        CCGAATG	CACATGTCT	Sample_28
        CCGAATG	TGTGCACGAT	Sample_29
        TTAGGCAG	CTAACGT	Sample_30
        TTAGGCAG	TCGGTACT	Sample_31
        TTAGGCAG	GATCGTTGT	Sample_32
        TTAGGCAG	AGCTACACTT	Sample_33
        TTAGGCAG	ACGCATT	Sample_34
        TTAGGCAG	GTATGCAT	Sample_35
        TTAGGCAG	CACATGTCT	Sample_36
        TTAGGCAG	TGTGCACGAT	Sample_37
        AACTCGTCG	CTAACGT	Sample_38
        AACTCGTCG	TCGGTACT	Sample_39
        AACTCGTCG	GATCGTTGT	Sample_40
        AACTCGTCG	AGCTACACTT	Sample_41
        AACTCGTCG	ACGCATT	Sample_42
        AACTCGTCG	GTATGCAT	Sample_43
        AACTCGTCG	CACATGTCT	Sample_44
        AACTCGTCG	TGTGCACGAT	Sample_45
        GGTCTACGTG	CTAACGT	Sample_46
        GGTCTACGTG	TCGGTACT	Sample_47
        GGTCTACGTG	GATCGTTGT	Sample_48
        GGTCTACGTG	AGCTACACTT	Sample_49
        GGTCTACGTG	ACGCATT	Sample_50
        GGTCTACGTG	GTATGCAT	Sample_51
        GGTCTACGTG	CACATGTCT	Sample_52
        GGTCTACGTG	TGTGCACGAT	Sample_53
        GATACCG	CTAACGT	Sample_54
        GATACCG	TCGGTACT	Sample_55
        GATACCG	GATCGTTGT	Sample_56
        GATACCG	AGCTACACTT	Sample_57
        GATACCG	ACGCATT	Sample_58
        GATACCG	GTATGCAT	Sample_59
        GATACCG	CACATGTCT	Sample_60
        GATACCG	TGTGCACGAT	Sample_61
        TCATGGTCAG	CTAACGT	Sample_62
        TCATGGTCAG	TCGGTACT	Sample_63
        TCATGGTCAG	GATCGTTGT	Sample_64
        TCATGGTCAG	AGCTACACTT	Sample_65
        TCATGGTCAG	ACGCATT	Sample_66
        TCATGGTCAG	GTATGCAT	Sample_67
        TCATGGTCAG	CACATGTCT	Sample_68
        TCATGGTCAG	TGTGCACGAT	Sample_69

#. The directory structure now looks something like this:

    .. code-block:: text

        .
        ├── AGTCGCTT+_R1_001.fastq.gz
        ├── AGTCGCTT+_R2_001.fastq.gz
        ├── CATTCGGT+_R1_001.fastq.gz
        ├── CATTCGGT+_R2_001.fastq.gz
        ├── CGAACTGT+
        │   └── CGAACTGT_internal_indexes.txt
        ├── CGAACTGT+_R1_001.fastq.gz
        ├── CGAACTGT+_R2_001.fastq.gz
        ├── demuxbyname.e631893
        ├── demuxbyname.o631893
        ├── demux.qsub
        ├── my_i7_indexes.txt
        ├── TCGGTTAC+_R1_001.fastq.gz
        ├── TCGGTTAC+_R2_001.fastq.gz
        ├── Undetermined_R1_001.fastq.gz
        └── Undetermined_R2_001.fastq.gz

#. Once that's done, we can start the demultiplexing with by creating a qsub script (``internal_dmux.qsub``) that looks like this:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <allocation>
        #PBS -l nodes=1:ppn=1
        #PBS -l walltime=12:00:00
        #PBS -q single
        #PBS -N radcap_demux_p1

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        process_radtags \
            -1 ../CGAACTGT+_R1_001.fastq.gz \
            -2 ../CGAACTGT+_R2_001.fastq.gz \
            -i gzfastq \
            -b CGAACTGT_internal_indexes.txt \
            --inline_inline \
            -o ./ \
            -c -q -r -t 140 -w 0.15 -s 10 \
            --renz_1 nheI \
            --renz_2 ecoRI \
            --adapter_mm 2 \
            --adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --retain-header

#. And, the directory structure looks like this:

    .. code-block:: text

        .
        ├── AGTCGCTT+_R1_001.fastq.gz
        ├── AGTCGCTT+_R2_001.fastq.gz
        ├── CATTCGGT+_R1_001.fastq.gz
        ├── CATTCGGT+_R2_001.fastq.gz
        ├── CGAACTGT+
        │   ├── CGAACTGT_internal_indexes.txt
        │   └── internal_dmux.qsub
        ├── CGAACTGT+_R1_001.fastq.gz
        ├── CGAACTGT+_R2_001.fastq.gz
        ├── demuxbyname.e631893
        ├── demuxbyname.o631893
        ├── demux.qsub
        ├── my_i7_indexes.txt
        ├── TCGGTTAC+_R1_001.fastq.gz
        ├── TCGGTTAC+_R2_001.fastq.gz
        ├── Undetermined_R1_001.fastq.gz
        └── Undetermined_R2_001.fastq.gz

#. Now, submit the job and let it run.  Proceed to do the same across the other plates of samples, adjusting the contents of the ``*_internal_indexes.txt`` file for whatever samples are in that plate.

    .. admonition:: Note

        It takes a while for stacks to write data to the files that are visible in the output directory you choose.  So, just be sure to give it some time to run (~10 minutes) before worrying too much about something wrong with how you set it up. The results also come somewhat slowly across all samples.

        Also be aware that the steps above usually take 3-6 hours to run.  If you have multiple plates of data, it would be sensible to setup jobs to demultiplex those, as well.

#. Once the data are demultiplexed, we need to remove duplicate reads.  Because we need to do this across many files (for each "plate"), we'll use GNU Parallel.  To get that process started, first make a new directory within each demultiplexed plate named ``process_radtags-removed`` and move all ``*.rem.1.fq.gz`` files there (alternatively, you can delete them).

    .. code-block:: bash

        mkdir removed
        mv *.rem.*.fq.gz ./process_radtags-removed/

#. You probably also want to check to see if any of the read files in the directory are **very** small (kB instead of MB).  You can do that with the following, which finds files < 25 MB in size. I would probably remove these samples from further consideration:

    .. code-block:: bash

        find . -maxdepth 1 -type f -size -25M

#. Next, make a new directory, ``duplicates-removed`` in the folder where you are working.  The directory structure will look like this:

    .. code-block:: text

        .
        ├── AGTCGCTT+_R1_001.fastq.gz
        ├── AGTCGCTT+_R2_001.fastq.gz
        ├── CATTCGGT+_R1_001.fastq.gz
        ├── CATTCGGT+_R2_001.fastq.gz
        ├── CGAACTGT+
        │   ├── CGAACTGT_internal_indexes.txt
        │   ├── duplicates-removed
        │   ├── internal_dmux.qsub
        │   ├── Sample_22.1.fq.gz
        │   ├── ...
        │   └── Sample_69.rem.2.fq.gz
        ├── CGAACTGT+_R1_001.fastq.gz
        ├── CGAACTGT+_R2_001.fastq.gz
        ├── demuxbyname.e631893
        ├── demuxbyname.o631893
        ├── demux.qsub
        ├── my_i7_indexes.txt
        ├── TCGGTTAC+_R1_001.fastq.gz
        ├── TCGGTTAC+_R2_001.fastq.gz
        ├── Undetermined_R1_001.fastq.gz
        └── Undetermined_R2_001.fastq.gz


#. Change into this new folder and generate an input file that will contain ``<sample_name>,<read1 path>,<read2 path>`` using the following:

    .. code-block:: bash

        for i in ../*.1.fq.gz; do b=`basename $i`; sample=${b%%.*}; echo "$sample,../$sample.1.fq.gz,../$sample.2.fq.gz" >> sample.list; done

#. We need to create a script that we will run with parallel that contains the code to remove duplicates.  Create a new file, ``clone_filter.sh`` and add to it the following:

    .. code-block:: bash

        #!/bin/bash

        SAMPL=$1
        READ1=$2
        READ2=$3

        # echo name of sample to stdout
        echo $SAMPL

        # remove PCR duplicates
        clone_filter \
            -P \
            -i gzfastq \
            --null-index \
            --oligo-len-2 8 \
            -1 $READ1 \
            -2 $READ2 \
            -D

#. Once that's done, you need to make it executable, so ``chmod +x clone_filter.sh``.
#. Now, create your qsub script in the same directory:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <allocation>
        #PBS -l nodes=1:ppn=1
        #PBS -l walltime=12:00:00
        #PBS -q single
        #PBS -N radcap_demux_p1

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # Number of Cores per job
        export CORES_PER_JOB=1

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
                -a sample.list \
                ./clone_filter.sh {$1} {$2} {$3}

