.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running RADcap Analysis
=======================

:Author: Brant C. Faircloth
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

#. To compile stacks when using LSU HPC, be sure to enable `module load gcc/6.4.0` in your `~/.modules` file 
#. Get stacks, configure (w/ home directory install), and install.  E.g.,

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

Steps
-----

#. Upload the relevant data to some location on @smic.  These should not have been demultiplexed in any way.  If you have multiple files (for some reason), you can combine the files together for READ1 and then combine the files together for READ2.

#. First, we need to demultiplex the data based on the i5 tags.  This may separate them by plate, or it may do something like separate them by row in a given plate.  To do this, you need to make a file of all the i5 index sequences that looks like this ``<tag-sequence>\t<plate/column name>``:

    .. code-block:: text

        CGAACTGT	109_2
        CATTCGGT	109_3
        TCGGTTAC	109_4
        AAGTCGAG	109_5
        TATCGGTC	109_6
        TATTCGCC	109_7
        GTATTGGC	109_8

#. Now, using that file, go ahead and demultiplex using ``process_radtags``.  The following command tells ``process_radtags`` that the data are paired (``-P``), the files are gzip, fastq file (``-i gzfastq``), that the barcodes (``-b``) are in ``i7_tags.txt``, to rescue (``-r``) barcodes that differ slightly, that the index is in the fastq header (``--index_null``), to not (yet) check for cut sites, and to retain header information:

    .. code-block:: bash

        process_radtags -P -i gzfastq \
            -1 ./Undetermined_S0_R1_001.fastq.gz \
            -2 ./Undetermined_S0_R2_001.fastq.gz \
	        -o ./ \
            -b ./i7_tags.txt \
            -r --index_null \
            --disable_rad_check --retain_header