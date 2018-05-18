.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

Random Computer Snippets
========================

Subsample reads for R1 and R2 using seqtk
-----------------------------------------

.. code-block:: bash
	
	for dir in /path/to/your/clearn/data/dir/from/illumiprocesser/*;
	do 
		RAND=$RANDOM;
		echo $RAND;
		for file in $dir/split-adapter-quality-trimmed/*-READ[1-2]*;
		do
			echo $file;
			seqtk sample -s $RAND 2000000 file | gzip > $file:t:r:r.SUBSAMPLE.fastq.gz
		done;
	done