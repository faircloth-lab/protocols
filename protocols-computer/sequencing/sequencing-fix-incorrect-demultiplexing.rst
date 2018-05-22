.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Fix Incorrect Demultiplexing
============================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.


Modification History
--------------------

See `Fix Incorrect Demultiplexing`_ 

.. _Fix Incorrect Demultiplexing: https://github.com/faircloth-lab/protocols/commits/master/protocols-computer/sequencing/sequencing-fix-incorrect-demultiplexing.rst

Purpose
-------

Sometimes, you get your data back from the sequencer, and you find that you have no data for some (many?) samples.  This usually results from the fact that you've made a mistake in getting the correct indexes to the sequencing facility.  Often, they are simply in the wrong orientation and you'll need to reverse complement one of the indexes to fix the problem and have your sequencing center demultiplex your data again.

Other times, the problem is a little more difficult.  This usually happens when a few indexes are mis-specified - so you get zero data for those samples, but lots of data from others.  This also usually results in *really* large ``Undetermined_R*_001.fastq.gz`` files.

Here's how I go about diagnosing the problem and potentially fixing it.


Steps
-----

1.  First, it's always good to get an idea of read counts for a given batch of samples.  If you have all of your `R1` and `R2` files in a directory, you can use something like the following to count reads in each file:

    .. code-block:: bash

        for i in *_R1_*; do echo $i; gunzip -c $i | wc -l; done

        009-03_Diprionidae_Neodiprion_edulicolus_R1_001.fastq.gz
        411108
        014-03B_Diprionidae_Neodiprion_sp__R1_001.fastq.gz
        784044
        016-03_Diprionidae_Neodiprion_ventralis_R1_001.fastq.gz
        1364944
        019-01_Diprionidae_Gilpinia_frutetorum_R1_001.fastq.gz
        278648
        044-03_Diprionidae_Neodiprion_autumnalis_R1_001.fastq.gz
        641604
        069-02_Diprionidae_Neodiprion_sertifer_R1_001.fastq.gz
        122256
        080-02_Diprionidae_Neodiprion_nanulus_nanulus_R1_001.fastq.gz
        271664
        090-03_Diprionidae_Neodiprion_nr._demoides_R1_001.fastq.gz
        354608

    These are **line** counts, so be sure to **divide these by 4 to get read counts**.  A pro-tip is that you can turn this into columns using the following find ``(.*)\n(.*)\n*`` and replace ``$1,$2\n`` commands work for your favorite text editor.


2.  You want to compare this list to what you expect, being aware of samples
that are either: (1) completely missing or (2) have very little data, like so:

    +----------------------------------------------------+------------+------------+
    | sample                                             | Line Count | Read Count |
    +====================================================+============+============+
    | myrmoborus_myotherinus_LSUMZ_5485_R1_001.fastq.gz  | 4          | 1          |
    +----------------------------------------------------+------------+------------+
    | myrmoborus_myotherinus_LSUMZ_74032_R1_001.fastq.gz | 8          | 2          |
    +----------------------------------------------------+------------+------------+
    | myrmoborus_myotherinus_LSUMZ_77634_R1_001.fastq.gz | 48         | 12         |
    +----------------------------------------------------+------------+------------+
    | myrmoborus_myotherinus_LSUMZ_907_R1_001.fastq.gz   | 48         | 12         |
    +----------------------------------------------------+------------+------------+
    | myrmoborus_myotherinus_MPEG_60007_R1_001.fastq.gz  | 4          | 1          |
    +----------------------------------------------------+------------+------------+

    These samples are likely some with incorrect indexes (we expected them to get lots of reads, but, in reality, they received very few).


3.  Take a peak into the undetermined file to get the sequencing machine name in
the header line:

    .. code-block:: bash

        gunzip -c Undetermined_R1_001.fastq.gz| less


    That looks like:

    .. code-block::text

        @J00138:114:HKMFKBBXX:7:1101:24951:1050 1:N:0:NCAAGACG+NATGTTCC
        ANCAAGTGATACATGTTCGATCTGTATGAATTCAGATAATTTTCTCATGTCGGGTAATATCTCACACTCAAGTATTTGCCAGTAACATTCTATCCGTCACTACATATGTTCCATATTTTATGTTCTTAGATCGGAAGAGCACACGTCTGAA
        +
        A#-FFJJJJJJJJ7FJFJJJJJJFFJJFJFFJJJJJJJJFFFJJAJJJJFJJJJJJJJJ<J<JJJJJFAFFJAJ-77<FFJJJJFJJJJJJJJJJJF7JJAJFFFJF7JJJJJFJFFFAFFJJJA-AAFJFJJ<FJF-<FFAJAFF7F)7<
        @J00138:114:HKMFKBBXX:7:1101:24992:1050 1:N:0:NATATCGA+NGATCTCG
        ANCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGCTCGCAAGGCTAACGACTCACACCACGACTACTAACAAAAATAATATATGCCTGAGTATGATACAACTAATAGCCGTCTTCATTTCCAAAC
        +
        A#AFFJJJJJFJJFJFAJJFFJFJJJJAFF-<FAJFF-FFJJFJJJJJJFAJFFJJF<7---FA<<7--7<AA--7--77AFF--7-7---7-7--A--7-7------777---77-A----7AA-7-7F<--7-7)-A7FFA--7A-7-7
        @J00138:114:HKMFKBBXX:7:1101:25093:1050 1:N:0:NGTGTCTG+NGATCTCG
        TNCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATAGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGC
        +
        A#AFFFJFJJJJJJJJFJJJJJJJFFJJFJJ<AAJ-FJJJF-FJAJJFA<FJJJJJJJJJJ7--FA<F--<<7JJFAF-A-<---FAJ7AA-77-77--<FJFA7<<-7F<7<7F7-A<FF<FFJ77<FFFJFJJ7<-7-<7<JJAAF---
        @J00138:114:HKMFKBBXX:7:1101:25134:1050 1:N:0:NAGAGCCA+NATACCAC
        CNTGTTAGCAGGTTGATTTGTGCACATTGGGACAATTAGTGGTTACTGTGAAACATTTTCCCTTCAGTGCTGTCAGTTCTGTATCTGCACTTTTCTTCCTTACCCTCACTGACATTGGCTCCCCTCCTTGGAGATCGGAAGAGCACACTTC
        +
        A#-AFJJJJJJ<A-FJJJJJJJJJJJJ7F7FJJJJJJJJAAAA-FJJJJJJJJFJJFJJJFJJJJJ<-<F7A<JFJ<FFJJJJJJJJJAJFJFJJAJAJJFFJJJJJFAJ7-AFJJ-7<J<-JA7F7FJ-)---AF7F7<A--AAAJF7-7

4. This was sequencer ``J00138``.  Now, parse out all the indexes in the
``Undetermined_R1_001.fastq.gz`` file and count them to see if you can see what happened.  Run the following:

    .. code-block:: bash

        gunzip -c  Undetermined_R1_001.fastq.gz | grep "^@J00138" | awk -F: '{print $NF}' | sort | uniq -c | sort -nr > R1_barcode_count.txt


5. Let's take a look in the ``R1_barcode_count.txt`` we just created.

    .. code-block:: bash

        less R1_barcode_count.txt

    Which looks like:

    .. code-block::text

        1833351 ACAGCTCA+TAGCGTCT
        1686592 ACAGCTCA+CATACCAC
        1537348 AAGAGCCA+CTACAGTG
        1342331 AAGAGCCA+AGCGTGTT
        1334012 AAGAGCCA+TTGCGAAG
        1217060 AAGAGCCA+TAGCGTCT
         988999 AAGAGCCA+TGGAGTTG
         972416 AAGAGCCA+GCTTCGAA
         940040 ACAGCTCA+CTACAGTG
         885663 AAAAAAAA+AGATCTCG <== this one is pretty useless
         768449 ACAGCTCA+TTGCGAAG
         716601 AAGAGCCA+ACCATCCA
         610761 AAGAGCCA+CATACCAC
         226926 AAATAAAA+AGATCTCG <== big decrease here


6. The first ~12 samples have a lot of reads associated with the given index sequences.  **You** now need to do a little detective work to see how things got screwed up and if these indexes match any/all of your missing samples.  Lots of times one of the indexes in the pair will be in the wrong orientation (so look at the revcomp of the index to help you solve the mystery).


7. Once you are pretty sure you have figured things out, download the tarball
for BBmap (https://sourceforge.net/projects/bbmap/).  Unzip that somewhere on
your machine.  This source has a really handy and fast script to parse out indexes, named ``demuxbyname.sh``.  After solving my missing sample mystery, I can parse those index combinations that I want to use into individual, index-specific ``R1`` and ``R2`` files using a command like the following.  The ``prefixmode=f`` command tells ``demuxbyname.sh`` to look at the suffix of the header line for the indexes specified by ``names=``:

    .. code-block:: bash

        ~/src/BBMap_37.33/demuxbyname.sh \
            prefixmode=f \
            in=../Undetermined_R1_001.fastq.gz \
            in2=../Undetermined_R2_001.fastq.gz \
            out=%_R1_001.fastq.gz \
            out2=%_R2_001.fastq.gz \
            names=AAGAGCCA+TTGCGAAG,AAGAGCCA+CATACCAC,AAGAGCCA+CTACAGTG,AAGAGCCA+TAGCGTCT,AAGAGCCA+TGGAGTTG,AAGAGCCA+AGCGTGTT,AAGAGCCA+ACCATCCA,AAGAGCCA+GCTTCGAA,ACAGCTCA+TTGCGAAG,ACAGCTCA+CATACCAC,ACAGCTCA+CTACAGTG,ACAGCTCA+TAGCGTCT

8. This will create a set of output files corresponding to ``R1`` and ``R2`` files for each of the index combinations.  On a pair of ~5 GB ``Undetermined_R*_001.fastq.gz`` files, this took about 250 seconds. That's fast. The output looks like:

    .. code-block:: bash

        -rw-r--r--. 1 bcf data  68207874 Jul  6 13:28 AAGAGCCA+ACCATCCA_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data  79802398 Jul  6 13:28 AAGAGCCA+ACCATCCA_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data 126637252 Jul  6 13:28 AAGAGCCA+AGCGTGTT_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 150255284 Jul  6 13:28 AAGAGCCA+AGCGTGTT_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data  57999953 Jul  6 13:28 AAGAGCCA+CATACCAC_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data  67958220 Jul  6 13:28 AAGAGCCA+CATACCAC_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data 145783313 Jul  6 13:28 AAGAGCCA+CTACAGTG_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 170994098 Jul  6 13:28 AAGAGCCA+CTACAGTG_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data  92062354 Jul  6 13:28 AAGAGCCA+GCTTCGAA_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 109146094 Jul  6 13:28 AAGAGCCA+GCTTCGAA_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data 114929942 Jul  6 13:28 AAGAGCCA+TAGCGTCT_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 136107891 Jul  6 13:28 AAGAGCCA+TAGCGTCT_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data  92648130 Jul  6 13:28 AAGAGCCA+TGGAGTTG_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 109343054 Jul  6 13:28 AAGAGCCA+TGGAGTTG_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data 125381989 Jul  6 13:28 AAGAGCCA+TTGCGAAG_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 148050908 Jul  6 13:28 AAGAGCCA+TTGCGAAG_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data 157282829 Jul  6 13:28 ACAGCTCA+CATACCAC_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 184159162 Jul  6 13:28 ACAGCTCA+CATACCAC_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data  88803030 Jul  6 13:28 ACAGCTCA+CTACAGTG_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 104882949 Jul  6 13:28 ACAGCTCA+CTACAGTG_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data 170820756 Jul  6 13:28 ACAGCTCA+TAGCGTCT_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data 200969341 Jul  6 13:28 ACAGCTCA+TAGCGTCT_R2_001.fastq.gz
        -rw-r--r--. 1 bcf data  72785440 Jul  6 13:28 ACAGCTCA+TTGCGAAG_R1_001.fastq.gz
        -rw-r--r--. 1 bcf data  86357340 Jul  6 13:28 ACAGCTCA+TTGCGAAG_R2_001.fastq.gz


   Each of the resulting files corresponds to the `R1` and `R2` reads for a given index combination.  There is also no error correction going on here (which is just fine by me).
