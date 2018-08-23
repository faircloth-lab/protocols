.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Demultiplexing a Sequencing Run
===============================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.


Modification History
--------------------

See `Demultiplexing a Sequencing Run`_ 

.. _Demultiplexing a Sequencing Run: https://github.com/faircloth-lab/protocols/commits/master/protocols-computer/sequencing/sequencing-demultiplex-a-run.rst


Purpose
-------

Often, we combine MANY libraries together in a single sequencing run (moreseo even now that NovaSeqs are online).  Once sequenced, the data generally need to be demultiplexed by their index sequences into something approximating the sample names that you want to associate the sequence data with.  You can generally do this one of two ways: (1) directly demultiplex to named files using `bcl2fastq` from Illumina or (2) You can receive `Undetermined` files from the sequencing center, and demultiplex those based on the index calls in the header of the sequence, for example:

    .. code-block:: text

        @A00484:41:H3G5VDRXX:1:1101:1018:1000 1:N:0:GGCGTTAT+CCTATTGG
                                                    ^^^^ indexes ^^^^

Steps
-----

1. Download sequence data from provider.  They will usually tell you how to do this - either with ``wget`` or with ``sftp``

2. If ``sftp`` and using a private certificate, you need to get the certificate info into a file (e.g. ``my.key``), then:
   
    .. code-block:: bash

        chmod 0600 my.key
        sftp -i /path/to/my.key user@sftp.some.edu

3. Things will take a long time to download.
4. Once downloaded, be sure to get/check md5sums of files against what provider gives you (often, these ``.md5`` files are part of the download).  These help you make sure that the downloads were not corrupted while downloading (which can happen with big files).

    .. code-block:: bash

        for i in *.md5; do md5sum -c $i; done

5. You may want to count the read numbers in the file.  You can check this number (times two) against the output below to make sure you're apples-to-apples on the overall count of reads. This can also take a long time.
   
    .. code-block:: bash

        gunzip -c Undetermined_S0_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}'

5. Get list of barcodes from the users who shared your run.  To make your life easy, they should provide a spreadsheet that's setup correctly, w/ both ``i5`` and ``i7`` names **and** sequences, as well as the forward **and** reverse complement of the ``i7`` sequences.

6. You're going to need to create combinations of barcode from the list provided by the users sharing a run. Before you do this, it's often easiest to peek inside the ``R1`` (or ``R2``) file from the run (typically named ``Undetermined_S0_L001_R1_001.fastq.gz``) to have a look at the index sequences reported and to make sure which of the forward or reverse complement of ``i7`` you need to use for a given platform.  To peek inside, use:

    .. code-block:: bash

        gunzip -c Undetermined_S0_L001_R1_001.fastq.gz | less

7. Then, enter "search mode" in ``less`` by typing ``/``. I then generally search for an ``i5`` index that will be prevalent in the run (if everything is equimolar, just pick one), then look to see what the ``i7`` sequence of the corresponding index is.  I compare that to my spreadsheet of forward and reverse indices for the ``i7`` position, and then I know what I need to do across all the ``i7`` indexes.

8. Once you're happy with that, create a file of barcodes, e.g. ``my_barcodes.txt`` where you have the indexes in the correct order, here ``reverse_comp(i7)+i5``:

    .. code-block:: text

        TTACCGAG+TCGTCTGA
        TTACCGAG+CATGTGTG
        TTACCGAG+TCTAGTCC
        TTACCGAG+AAGGCTCT
        TTACCGAG+AACCAGAG
        TTACCGAG+ACTATCGC
        TTACCGAG+GTCCTAAG
        TTACCGAG+TGACCGTT
        GTCCTAAG+TCGTCTGA
        GTCCTAAG+CATGTGTG
        GTCCTAAG+TCTAGTCC
        GTCCTAAG+AAGGCTCT
        GTCCTAAG+AACCAGAG
        GTCCTAAG+ACTATCGC
        GTCCTAAG+GTCCTAAG
        GTCCTAAG+TGACCGTT

9. You can use that file and ``demuxbyname.sh`` from BBMap (here v38.22) to demultiplex paired files of ``Unknown`` reads into resulting files that will be labelled with their respective indexes:

    .. code-block:: bash

        ~/src/BBMap_38.22/demuxbyname.sh \
            prefixmode=f \
            in=../Undetermined_S0_L001_R1_001.fastq.gz \
            in2=../Undetermined_S0_L001_R2_001.fastq.gz \
            out=%_R1_001.fastq.gz \
            out2=%_R2_001.fastq.gz \
            outu=Undetermined_R1_001.fastq.gz \
            outu2=Undetermined_R2_001.fastq.gz \
            names=../index_sequences.txt

10. For a NovaSeq S1 run, this took about 2.3 hours and produced, as output:
    
    .. code-block:: bash

        Input is being processed as paired
        Time:               8250.056 seconds.
        Reads Processed:    2196226240          266.21k reads/sec
        Bases Processed:    331630162240        40.20m bases/sec
        Reads Out:    4080560900
        Bases Out:    616164695900


11. Once you have these files, you are almost there.  The easiest thing to do to get all the files renamed is to create another, tab-delimited list thas has the index combinations in column 1 and the names you want for the file in column 2.  Name this file ``temp-names.txt``.  This file looks something like:

    .. code-block:: bash

        TTACCGAG+TCGTCTGA       Molothrus_ater_LSUMZ441_NT
        TTACCGAG+CATGTGTG       Molothrus_ater_Tulane2906_NT
        TTACCGAG+TCTAGTCC       Molothrus_ater_LSUMZ2237_NT
        TTACCGAG+AAGGCTCT       Molothrus_ater_LSUMZ2241_NT
        TTACCGAG+AACCAGAG       Molothrus_ater_LSUMZ5504_NT
        TTACCGAG+ACTATCGC       Molothrus_ater_LSUMZ13890_NT 
    

12. Now you can rename the files in your set by running:
    
    .. code-block:: bash

        while IFS=$'\t' read -r column1 column2; do
            mv ${column1}_R1_001.fastq.gz ${column2}_${column1}_R1_001.fastq.gz;
            mv ${column1}_R2_001.fastq.gz ${column2}_${column1}_R2_001.fastq.gz;
        done < "temp-names.txt"

13. This keeps the original tag name in the name of the raw data, but prepends the name you want from the 2nd column of the tab-delimited file, above.  The files names now look something like:
    
   .. code-block:: bash

    Molothrus_ater_LSUMZ13890_NT_TTACCGAG+ACTATCGC_R1_001.fastq.gz
    Molothrus_ater_LSUMZ13890_NT_TTACCGAG+ACTATCGC_R2_001.fastq.gz
    Molothrus_ater_LSUMZ160263_NT_GTCCTAAG+GTCCTAAG_R1_001.fastq.gz
    Molothrus_ater_LSUMZ160263_NT_GTCCTAAG+GTCCTAAG_R2_001.fastq.gz
    Molothrus_ater_LSUMZ160263_P_GAAGTACC+TGACCGTT_R1_001.fastq.gz
    Molothrus_ater_LSUMZ160263_P_GAAGTACC+TGACCGTT_R2_001.fastq.gz
    Molothrus_ater_LSUMZ160263_PW_CAGGTATC+AACCAGAG_R1_001.fastq.gz
    Molothrus_ater_LSUMZ160263_PW_CAGGTATC+AACCAGAG_R2_001.fastq.gz
    Molothrus_ater_LSUMZ160264_NT_GTCCTAAG+TGACCGTT_R1_001.fastq.gz
    Molothrus_ater_LSUMZ160264_NT_GTCCTAAG+TGACCGTT_R2_001.fastq.gz 
