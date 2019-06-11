.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _AssemblyWithCanu:

Assembly With Canu
==================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Assembly With Canu`_.

.. _Assembly With Canu: https://github.com/faircloth-lab/protocols/commits/master/protocols-computer/assembly/assembly-with-canu.rst


Purpose
-------

There are several options to assemble PacBio long-read data, and one of those (potentially the easier to install) is canu_ (another is to the the `SMRTAnalysis <https://www.pacb.com/support/software-downloads/>`_  pipeline and/or `pb-assembly <https://github.com/PacificBiosciences/pb-assembly>`_ ).  canu_ works reasonably well on @QB2 - I've just generally learned that it's easier to run in single-threaded mode rather than try to make the grid mode work (it seems as if grid mode most likely will NOT work on the queueing system that we use.

Steps
-----

#. Because canu_ is compute intensive, the following steps have been documented assuming you are using @QB2
#. Compile canu_ according to :ref:`CompilingCanu`
#. Create a `pacbio` environment for conda (after installing miniconda_ and configuring for bioconda_)

    .. code-block:: bash

        conda create -n pacbio python=2.7 bam2fastx
    

#. Create a working directory for your data (here, I'm just using the _Arabidopsis_ test data):

    .. code-block:: text

      mkdir arabidopsis-pacbio && cd arabidopsis-pacbio

#. Download *Arabidopsis* test data from PacBio.  Be sure to get the ``*.pbi`` files because we need them to convert the ``bam`` data to ``fastq`` format

    .. code-block:: text

      wget -P pacbio-raw https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/1_A01_customer/m54113_160913_184949.subreads.bam
      wget -P pacbio-raw https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/1_A01_customer/m54113_160913_184949.subreads.bam.pbi

      wget -P pacbio-raw https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/3_C01_customer/m54113_160914_092411.subreads.bam
      wget -P pacbio-raw https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/3_C01_customer/m54113_160914_092411.subreads.bam.pbi


#. canu_ requires data in ``fastq`` format, so convert each ``bam`` file to ``fastq``.

    .. code-block:: bash

      #!/bin/bash
      #PBS -q single
      #PBS -A <allocation>
      #PBS -l walltime=06:00:00
      #PBS -l nodes=1:ppn=2
      #PBS -V
      #PBS -N bam_to_fastq
      #PBS -o bam_to_fastq.out
      #PBS -e bam_to_fastq.err

      # load the parallel module to run files in parallel (up to 4 cores in single queue)
      module load gnuparallel/20170122

      # activate our conda env
      source activate pacbio


      cd $PBS_O_WORKDIR
      mkdir pacbio-fastq && cd pacbio-fastq
      find ../pacbio-raw/ -name *.bam | parallel "bam2fastq -o {/.} {}"

   .. note:: You may need to adjust queues and cores to suit your needs.  Here, I'm using the ``single`` queue because I only have 2 files to convert and we can use up to 4 CPUs in ``single``.  Also note that you may need to adjust the time needed for each run - particularly for larger bam files you are converting.


#. Once those data are converted, we can kick off the canu_ assembly job.  Again, I've found that we need to keep these assembly jobs "local", meaning that we're not going to run in grid mode.  However, you do want to run them using the ``bigmem`` queue on @QB2.  Also note here that we're redirecting `stdout` and `stderr` to files - we're doing this so that we can check on job status as the runs go along (since the queuing system typically keeps these in temp files until the end of the run):

    .. code-block:: bash

      #!/bin/bash
      #PBS -q bigmem
      #PBS -A <allocation>
      #PBS -l walltime=72:00:00
      #PBS -l nodes=1:ppn=48
      #PBS -V
      #PBS -N canu_config
      #PBS -o canu_config.out
      #PBS -e canu_config.err
      #PBS -m abe
      #PBS -M brant@faircloth-lab.org

      module load gcc/6.4.0
      module load java/1.8.0

      cd $PBS_O_WORKDIR
      mkdir -p canu-assembly && cd canu-assembly

      canu \
          -p arabidopsis \
          -d arabidopsis-pacbio \
          genomeSize=123m \
          useGrid=false \
          -pacbio-raw ../pacbio-fastq/*.fastq.gz 1>canu-assembly.stdout 2>canu-assembly.stderr
   
   .. warning:: You will need to adjust the genome size of your organism in the above to something that's suitable.  Gb-size genome size is set using ``genomeSize=1.1g``, which would be appropriate for a bird.

#. If you need to restart the job at any time (e.g., you run out of walltime, which is likely), you may want to rename ``canu-assembly.stdout`` and ``canu-assembly.stderr``, so they are not overwritten:

    .. code-block:: bash
    
      for i in canu-assembly/*.std*; do mv $i $i.old; done

#. Then you simply need to resubmit the ``qsub`` script and the job will restart from where it last started. 