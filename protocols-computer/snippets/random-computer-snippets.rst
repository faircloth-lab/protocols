.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

Random Computer Snippets
========================

All of the following assume that you are using the Z shell (zsh).  These may or may not work in BASH.

Subsample reads for R1 and R2 using seqtk
-----------------------------------------

.. code-block:: bash
	
	READS=2000000
	for dir in /path/to/your/clearn/data/dir/from/illumiprocesser/*;
	do 
		RAND=$RANDOM;
		echo $RAND;
		for file in $dir/split-adapter-quality-trimmed/*-READ[1-2]*;
		do
			echo $file;
			seqtk sample -s $RAND $file $READS | gzip > $file:t:r:r.SUBSAMPLE.fastq.gz
		done;
	done

Download data for multiple files from NCBI SRA
----------------------------------------------

First, create a list of SRRs in a file, `sra-records.txt`, that looks something like:

.. code-block:: text

	SRR453553
	SRR453556
	SRR453559
	SRR453277
	SRR453409
	SRR453550
	SRR452995
	SRR453269
	SRR453270
	SRR453274
	SRR453263

`fastq-dump` is slow AF over the internet. It's faster to download the SRA file using wget, then parse that into fastq:

.. code-block:: bash

	for record in `cat sra-records.txt`; 
	do
		echo $record;
		wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${record[1,6]}/$record/$record.sra;
		fastq-dump -I --split-files --gzip $record.sra;
		rm $record.sra;
	done

This is still pretty darn slow.  We can parallelize using GNU `parallel`.  After creating `sra-records.txt` to hold the list of files that we need, create a script named `get_sra.sh` with the following contents:

.. code-block:: bash

	#!/bin/zsh

	echo "starting $1";
	wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${1[1,6]}/$1/$1.sra;
	fastq-dump -I --split-files --gzip $1.sra;
	rm $1.sra;
	echo "finished $1";

Once, that's completed, make it executable (`chmod 0755 get_sra.sh`). Now, you can run this on 12 cores (`--jobs 12`) using:

.. code-block:: bash

	parallel --progress --joblog logfile --jobs 12 -a sra-records.txt ./get_sra.sh {$1}

This will download the files from SRA and convert them to fastq.gz in parallel.

