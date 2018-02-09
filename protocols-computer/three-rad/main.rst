.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


3RAD Demultiplexing and Analysis
================================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

.. program-output:: git log --graph --oneline --decorate -- protocols-computer/three-rad/main.rst

Purpose
-------


Preliminary Steps
-----------------


#. To compile stacks when using LSU HPC, be sure to enable +gcc-4.9.2 in your `~/.soft` file 
#. Get stacks, configure (w/ home directory install), and install.  E.g.,

    .. code-block:: bash

        wget http://catchenlab.life.illinois.edu/stacks/source/stacks-1.48.tar.gz
        tar -czvf stacks-1.48.tar.gz

        ./configure --prefix=/home/brant/
        make
        make install

Steps
-----

#. Upload relevant raw-read data to supermike.  These should have been demultiplexed by OUTER i5 and i7 tags already (this usually means they are already in plates).  We'll put these in ``$HOME/threerad/raw-reads``.
#. You may want to have each plate's-worth of data in it's own folder, then work within each folder across all plates.
#. Within each plate of samples, here are the internal tags (for the NheI + EcoRI combination):
   
    .. code-block:: text

        iTru_NheI_R1_stub_A A   CCGAAT
        iTru_NheI_R1_stub_B B   TTAGGCA
        iTru_NheI_R1_stub_C C   AACTCGTC
        iTru_NheI_R1_stub_D D   GGTCTACGT
        iTru_NheI_R1_stub_E E   GATACC
        iTru_NheI_R1_stub_F F   AGCGTTG
        iTru_NheI_R1_stub_G G   CTGCAACT
        iTru_NheI_R1_stub_H H   TCATGGTCA
                
        iTru_EcoRI_R2_1 1   CTAACGT
        iTru_EcoRI_R2_2 2   TCGGTACT
        iTru_EcoRI_R2_3 3   GATCGTTGT
        iTru_EcoRI_R2_4 4   AGCTACACTT
        iTru_EcoRI_R2_5 5   ACGCATT
        iTru_EcoRI_R2_6 6   GTATGCAT
        iTru_EcoRI_R2_7 7   CACATGTCT
        iTru_EcoRI_R2_8 8   TGTGCACGAT
        iTru_EcoRI_R2_9 9   GCATCAT
        iTru_EcoRI_R2_10    10  ATGCTGTT
        iTru_EcoRI_R2_11    11  CATGACCTT
        iTru_EcoRI_R2_12    12  TGCAGTGAGT


#. Before proceeding, you need to make sure that the barcode file that you are using accounts for the fact that the "left side" index needs a "G" added to the 5' end of the index sequence and the "right side" index sequence needs a "T" added to the 5' end of the index sequence.  In this way, you get:
   
    .. code-block:: text

        [Left Side]

        CCGAAT      CCGAATG
        TTAGGCA     TTAGGCAG
        AACTCGTC    AACTCGTCG
        GGTCTACGT   GGTCTACGTG

        [Right Side]

        CTAACG      CTAACGT
        TCGGTAC     TCGGTACT
        GATCGTTG    GATCGTTGT
        AGCTACACT   AGCTACACTT

#. For the NheI + EcoRI combination, I've done all this in a spreadsheet with locked cells (stacks-worksheet.xlsx), so that you only need to enter your sample names for each well.  You can download that spreadsheet from https://www.dropbox.com/s/lr6p7k5894wghux/stacks-worksheet.xlsx?dl=0.

#. Enter your sample name details in the appropriate column of the 2nd worksheet tab.

#. That results produces a stacks barcode file, which looks like (excerpt of a few lines):

    .. code-block:: text

        CCGAATG CTAACGT Sample1
        CCGAATG GATCGTTGT       Sample2
        CCGAATG ACGCATT Sample3
        CCGAATG CACATGTCT       Sample4
        CCGAATG GCATCAT Sample5
        CCGAATG CATGACCTT       Sample6
        TTAGGCAG        CTAACGT Sample7
        TTAGGCAG        GATCGTTGT       Sample8


#. Before proceeding, you need to manually create an output directory:

    .. code-block:: bash
        
        mkdir test-out

#. Now, you can demultiplex samples using something like the following (use -D if you want to keep discarded reads).  Note that the following assumed that ``process_radtags`` is in your ``$PATH``.  We are running this from ``$HOME/threerad``:
    
    .. code-block:: text

        process_radtags -1 $HOME/threerad/raw-reads/Brant_1A_S55_R1_001.fastq.gz -2 $HOME/threerad/raw-reads/Brant_1A_S55_R2_001.fastq.gz \
            -i gzfastq \
            -b test-tags.txt --inline_inline \
            -o demultiplex \
            -c -q -r -t 140 -w 0.15 -s 10\
            --renz_1 xbaI \
            --renz_2 ecoRI \
            --adapter_mm 2 \
            --adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

#. This will run for some time and finally produce some output that looks like:

    .. code-block:: text

        2000000 total sequences
        6325    reads contained adapter sequence (0.3%)
        51868   ambiguous barcode drops (2.6%)
        50133   low quality read drops (2.5%)
        22074   ambiguous RAD-Tag drops (1.1%)
        1869600 retained reads (93.5%)

#.  The above should also place the data in ``$HOME/threerad/demultiplex``. After the file are demultiplexed, we have two choices to make - *De novo Analysis* or *Reference-based Analysis*

Reference-based Analysis
^^^^^^^^^^^^^^^^^^^^^^^^

#. On supermike, we should install a recent version of ``bwa``.  You can do that by building from source:
    
    .. code-block:: bash

        mkdir $HOME/src
        cd src
        wget $(wget -O - https://api.github.com/repos/lh3/bwa/releases | grep tarball_url| head -n 1 | cut -d '"' -f 4)
        tar -xzvf v0.*
        cd v0.*
        make

#. Make sure the binaries produced are in your ``$PATH``.

#. We also need a reasonably recent version of samtools, but you can use the version that is on Supermike for the following. If you build ``samtools``, make sure the binaries produced are in your ``$PATH``.

#. Now, we need to get our reference genome from wherever it is located, and upload it to supermike (or similar). As before, I will assume you are running everything in ``$HOME/work/threerad``.  Upload your genome to ``$HOME/work/threerad/genome/my_genome.fasta``

#. We need to index the reference genome that we are using, prior to running alignment (this may take some time):

    .. code-block:: bash

        cd $HOME/work/threerad/genome/
        bwa index my_genome.fasta

#. Go back to the top level of the directory where we are working:
   
    .. code-block:: bash
   
        cd $HOME/work/threerad/

#. Now, we have a couple of options.  If you have them, you can take advantage of multi-node HPC resources using ``GNU parallel`` to spread jobs across multiple nodes. Or, you can run everything on a single node (with multiple cores).  Generally, you want to take advantage of the HPC resources, so we'll cover that path first.


Running stacks on multiple HPC nodes
""""""""""""""""""""""""""""""""""""

#. We are going to be running stacks using ``GNU parallel`` (hereafter, ``parallel``).  This is installed on SuperMike and should be enabled by adding ``+gnuparallel-20161022-gcc-4.4.6`` to your ``~/.soft`` file.

#. The setup for using ``parallel`` generally requires: (1) a QSUB script to start the job, (2) a script that QSUB calls to run individual parts of each job, and (3) a list of files to be input from #1 to #2.

#. First we need to create a shell script that will run ``bwa`` for each of our samples.  We will do this in the same directory where we want the output to go.  So, first, create this directory using ``mkdir bwa-alignments``, then change into that directory.  Next, create a bash script in ``bwa-alignments`` named ``multi_bwa.sh`` that looks like:

    .. code-block:: bash
           
        #!/bin/bash

        ## set this manually
        CORES_PER_JOB=4

        ## DO NOT EDIT BELOW THIS LINE - this comes as input from GNU parallel on STDIN, via the sample.list file

        SAMPL=$(basename $1)
        READ1=$2
        READ2=$3
        GENOME=$4
        GENOME_NAME=$(basename $GENOME)

        ## Here are the specific commands we are running

        # run bwa and output BAM
        bwa mem -t $CORES_PER_JOB $GENOME $READ1 $READ2 | samtools view -bS - > $SAMPL.bam

        # filter BAM for duplicate reads, imperfect matches, and >5 SNPs per read (NM:i:[0-5], below)
        samtools view -h -q 25 -F 4 -F 256 $SAMPL.bam | grep -v XA:Z | grep -v SA:Z | awk '{if($0 ~ /^@/ || $6 ~ /140M/) {print $0}}' | grep -E '^@|NM:i:[0-5]\s' | samtools view -bS - > $SAMPL.q30.unique.perfect.bam

        # remove the unfiltered BAM file
        rm $SAMPL.bam
    
#. You can edit ``CORES_PER_JOB`` if you would like to allocate more cores to each ``bwa`` job.  After creating this file, we need to make it executable by running ``chmod 0755 multi_bwa.sh``.
#. Next, create the QSUB script for the ``multi_bwa.sh`` script named ``bwa_run.qsub``.  It should look like the following (be sure that you replace ``<allocation_name>`` with your allocation name and update ``CORES_PER_JOB`` if you changed that, above).  Also, as written, below, this requests **12** nodes of **16** cores each, for a total of **192** cores.  We'll use 4 cores per job, so this will run **192/4 = 48** samples simultaneously.  Be sure to adjust the number of nodes requested if you have substantially more or fewer samples than this:

    .. code-block:: bash
           
        #PBS -A <allocation_name>
        #PBS -l nodes=12:ppn=16
        #PBS -l walltime=2:00:00
        #PBS -q checkpt
        #PBS -N multi_bwa

        # Number of Cores per job (needs to be multiple of 2)
        export CORES_PER_JOB=4

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
                ./multi_bwa.sh {$1} {$2} {$3} {$4}

#. Finally, you need to create the ``sample.list`` that will be read by ``bwa_run.qsub`` and passed to ``multi_bwa.sh``.  The easiest way to do this is to: (1) note the path to your ``bwa``-indexed genome, then (2) run (be sure to use the ``$PATH`` to *your* genome):
   
    .. code-block:: bash
   
        cd $HOME/work/threerad
        GENOME=$HOME/work/threerad/genome/my_genome.fasta
        ls -d $PWD/raw-reads/*.1.fq.gz | sed -E "s/(.*).1.fq.gz/\1,\1.1.fq.gz,\1.2.fq.gz/" | sed -E "s|.*|&,$GENOME|" > bwa-alignments/sample.list

#. This will create the file ``bwa-alignments/sample.list``, which will contain the following columns of information:

    .. code-block:: text
       
        1. path to your sample name
        2. path to your Read 1 files
        3. path to your Read 2 files
        4. path to your indexed genome

#. After that's all done, you can submit the QSUB script by running ``qsub bwa_run.qsub``.  You jobs should start once the queue has room, and they should not take too long to run (the 2 hour queue-time is sufficient to align several hundred MBs of data for each sample).
   
#. Next, we need to run ``pstacks`` against all of the BAM files we just created.  to do that, first ``cd $HOME/work/threerad``, then ``mkdir stacks`` and ``cd stacks``.  In this directory, we need to create a shell script to run pstacks ``multi_pstacks.sh``, a QSUB file to submit that job (``pstacks_run.qsub``), and a ``sample-bam.list``.  Create ``multi_pstacks.sh`` first so that it looks like (again, you can change ``CORES_PER_JOB``, but be sure to also change this in ``pstacks_run.qsub``, if you do):

    .. code-block:: bash

        #!/bin/bash

        ## set this manually
        CORES_PER_JOB=4

        ## DO NOT EDIT - this comes as input from GNU parallel on STDIN, via the sample.list file

        INTEGER=$1
        BAM=$2

        ## Here are the specific commands we are running
        # run pstacks
        pstacks -p $CORES_PER_JOB -t bam -m 3 -i $INTEGER -f $BAM -o ./

#. Make this executable by running ``chmod 0755 multi_pstacks.sh``.  Now, create a ``pstacks_run.qsub`` that looks like:
   
    .. code-block:: bash

        #PBS -A <allocation_name>
        #PBS -l nodes=3:ppn=16
        #PBS -l walltime=2:00:00
        #PBS -q checkpt
        #PBS -N multi_pstacks

        # Number of Cores per job (needs to be multiple of 2)
        export CORES_PER_JOB=4

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
                -a sample-bam.list \
                ./multi_pstacks.sh {$1} {$2}

#. Note that the above uses fewer nodes (the jobs are less compute intense).  Finally, we need to create ``sample-bam.list``, which we can do by running:

    .. code-block:: bash
       
        cd $HOME/work/threerad
        ls $PWD/bwa-alignments/*.bam | awk '{printf "%d,%s\n", NR, $0}' > stacks/sample-bam.list

#. This will create the file ``stacks/sample-bam.list``, which will contain the following columns of information:

    .. code-block:: text
   
        1. an integer value, unique to each sample
        2. the path to each sample's BAM file, created above

#. Note that the first value of ``stacks/sample-bam.list`` needs to be an integer value unique to each sample.  Once all that is done, you can:
   
    .. code-block:: bash

        cd $HOME/work/threerad/pstacks
        qsub pstacks_run.qsub

#. Now, we need to run ``cstacks`` against the resulting data.  We can't spread this job across multiple nodes (but we can use multi-threading).  First, ``cd $HOME/work/stacks``.  Then create a ``cstacks.sh`` script to run ``cstacks`` that looks like the following:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <your allocation>
        #PBS -l nodes=1:ppn=16
        #PBS -l walltime=12:00:00
        #PBS -q checkpt
        #PBS -N multiple_bwa

        #move into the directory containing this script
        cd $PBS_O_WORKDIR

        STACKS=$HOME/work/threerad/stacks

        # Create a list of file to supply to cstacks

        samples=""
        for file in $STACKS/*.tags.tsv.gz;
        do 
            prefix=$(echo $file | sed -E "s/.tags.tsv.gz//");
            samples+="-s $prefix ";
        done

        # Build the catalog; the "&>>" will capture all output and append it to the Log file.

        cstacks -g -p 16 -b 1 -n 1 -o ./ $samples &>> ./cstacks.log

#. Now that we've run ``cstacks``, we need to run ``sstacks``.  We can go back to running this in parallel, as before.  In ``$HOME/work/stacks``, create a script to run ``sstacks`` named ``multi_sstacks.sh``:

    .. code-block:: bash

        #!/bin/bash

        ## set this manually
        CORES_PER_JOB=4

        ## DO NOT EDIT - this comes as input from GNU parallel on STDIN, via the sample.list file

        SAMPLE=$1

        # run sstacks
        sstacks -g -p 4 -b 1 -c ./batch_1 -s $SAMPLE -o ./ &>> sstacks.log

#. Make the above script executable with ``chmod 0755 multi_sstacks.sh``. Next, create a QSUB file named ``sstacks_run.qsub``:
   
    .. code-block:: bash

        #PBS -A <allocation_name>
        #PBS -l nodes=3:ppn=16
        #PBS -l walltime=2:00:00
        #PBS -q checkpt
        #PBS -N multi_pstacks

        # Number of Cores per job (needs to be multiple of 2)
        export CORES_PER_JOB=4

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --progress \
                --joblog logfile.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a samples.list \
                ./multi_sstacks.sh {$1}

#. As before, you can change the number of nodes and number of ``CORES_PER_JOB``, but you need to do that in **both** files.  Finally, assuming we are in ``$HOME/work/stacks``, we need to create a ``samples.list`` file containing a list of all the samples we want to process:
   
    .. code-block:: bash

        ls $PWD/*.tags.tsv.gz | sed -E "s/.tags.tsv.gz//" | grep -v "catalog" > samples.list


#. Finally, we can run the job using ``qsub sstacks_run.qsub``.

#. Once ``sstacks`` has finished running, we can run ``populations``.  The following script runs populations in a manner identical to what is run during the final stages of the ``ref_map.pl`` script that is distributed with stacks:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <your allocation>
        #PBS -l nodes=1:ppn=16
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N stacks_pop_std

        #move into the directory containing this script
        cd $PBS_O_WORKDIR

        STACKS=$HOME/work/threerad/stacks

        # Run the populations code
        echo "=========`date`=========" >> populations.log
        populations -b 1 -P $STACKS -s -t 16 &>> populations.log

#. In similar fashion, we can create a VCF file from all of the data with the following script (``population_vcf.sh``).  Here, we're specifying a very low ``--min_maf`` to ``populations``, so that stacks will compute some allele frequency stats for each variant (we can filter this later with ``vcftools`` or similar):
   
    .. code-block:: bash

        #!/bin/bash
        #PBS -A <your allocation>
        #PBS -l nodes=1:ppn=16
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N stacks_pop_vcf

        #move into the directory containing this script
        cd $PBS_O_WORKDIR

        STACKS=$HOME/work/threerad/stacks

        # Run the populations code
        echo "=========`date`=========" >> populations.log
        populations -b 1 -P $STACKS --min_maf 0.02 --vcf -t 12 &>> populations.log


.. Running stacks on a single node
.. """""""""""""""""""""""""""""""



Additional Filtering
^^^^^^^^^^^^^^^^^^^^

#. We can filter the data in a number of ways once we have a VCF file:
   
    .. code-block:: bash

        # get a count of the missingness, by individual - outputs results into 'out.imiss'
        vcftools --vcf batch_1.vcf --missing-indv

        # filter the vcf file for missing data (this is no missing data for any individual).  --recode-INFO-all keeps all INFO fields
        vcftools --vcf batch_1.vcf --max-missing 1 --recode --recode-INFO-all --out no_missing

        # do same as above, but allow some missingness and filter for MAF ≥ 0.1
        vcftools --vcf batch_1.MAF.vcf --max-missing 0.95 --recode --maf 0.1 --recode-INFO-all --out 5p_missing_0.2-maf


