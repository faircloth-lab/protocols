.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running GATK in Parallel
========================

:Author: Jessie Salter, Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running GATK in Parallel`_ 

.. _Running GATK in Parallel: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-gatk-parallel.rst

Purpose
-------

Nothing here yet.


Preliminary Steps
-----------------

#. Install miniconda following the `instructions for bioconda <https://bioconda.github.io/user/install.html>`_.

#. Download a version of GATK from `their website <https://gatk.broadinstitute.org/hc/en-us/>`_.

#. Unzip that package.  Ensure that within the package there is a ``gatkcondaenv.yml`` file.

#. Now, we'll build a conda environment from the ``yml`` file.  This will take a little while, so you probably want to start an interactive job on the HPC, so your job doesn't die due to time limits:

    .. code-block:: bash

        qsub -I -l walltime=02:00:00,nodes=1:ppn=20 -A <allocation>

        # once the interactive session starts, navigate to the GATK
        # package and install
        cd <location of the unziped gatk package>
        conda env create -n gatk -f gatkcondaenv.yml

#. Link in the ``gatk`` binary:

    .. code-block:: bash

        conda activate gatk
        # find out where python lives
        which python
        # change to the bin directory in this environment and
        # link to the GATK wrapper which is back in the GATK package
        ln -s <path to gatk wrapper>

#. Theoretically, you will also want ``bwa`` and ``samtools`` in this directory to make your life easier.  You can install those with:

    .. code-block:: bash

        conda activate gatk
        conda install bwa samtools=1.9

#. That said, the above can be extremely slow. You may want to want to create another environment with an up-to-date ``bwa`` and ``samtools``.  This is usually much faster.  Be aware of the approach you take because it is important, later, in terms of how you call ``GATK`` relative to ``bwa`` or ``samtools``.

    .. code-block:: bash

        conda create -n mapping bwa samtools=1.9

#. Once that's done, we probably also want to go trimmomatic, so we can perform read trimming.  You can do that using ``wget`` to download the file and putting the trimmomatic jar file somewhere (``$HOME/jar/Trimmomatic-0.39/trimmomatic-0.39.jar``):

    .. code-block:: bash

        mkdir $HOME/jar
        cd $HOME/jar
        wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
        unzip Trimmomatic-0.39.zip

Steps
-----

#. When you get started, your read files are going to look something like this. Note that there is a similar pattern to the read files where the first part of each read name is the same up to the ``*.1.fastq.gz`` and ``*.2.fastq.gz`` - this is pretty common and we are going to take advantage of that:

    .. code-block:: bash

        ./
        └── raw-reads
            ├── HC2HMDSXX.1.AGCTTT.unmapped.1.fastq.gz
            ├── HC2HMDSXX.1.AGCTTT.unmapped.2.fastq.gz
            ├── HC2HMDSXX.1.AGGAAT.unmapped.1.fastq.gz
            ├── HC2HMDSXX.1.AGGAAT.unmapped.2.fastq.gz
            ├── HC2HMDSXX.1.AGTGCC.unmapped.1.fastq.gz
            ├── HC2HMDSXX.1.AGTGCC.unmapped.2.fastq.gz
            ├── HC2HMDSXX.1.AGTTCC.unmapped.1.fastq.gz
            ├── HC2HMDSXX.1.AGTTCC.unmapped.2.fastq.gz
            ├── HC2HMDSXX.1.ATCCGC.unmapped.1.fastq.gz
            ├── HC2HMDSXX.1.ATCCGC.unmapped.2.fastq.gz
            ├── HC2HMDSXX.1.ATGACT.unmapped.1.fastq.gz
            └── HC2HMDSXX.1.ATGACT.unmapped.2.fastq.gz

#. We need to generate an input file of the common portion of the read name, plus the count of threads that we want to use to trim each set of reads (you'll need to determine this, here I'll use 4 threads per set of read file).  You can do this in any number of ways, including using an external editor, excel, sed, a bash loop, etc.  Here's a ``sed`` example:

    .. code-block:: bash

        ls raw-reads/*.1.fastq.gz | sed -r "s/(.*).1.fastq.gz/\1,4/" > files-to-trim.txt

#. Now, we will make a script to trim these files, ``trimmomatic-sub.sh``

    .. code-block:: bash

        #!/bin/bash

        # name of output folder
        OUTPUT=clean-reads

        # DONT EDIT BELOW
        READ1=$1.1.fastq.gz
        READ2=$1.2.fastq.gz
        PREFIX=`basename $1`
        CLEAN_PAIRED_READ1=$OUTPUT/$PREFIX.1.clean.paired.fastq.gz
        CLEAN_PAIRED_READ2=$OUTPUT/$PREFIX.2.clean.paired.fastq.gz
        CLEAN_UNPAIRED_READ1=$OUTPUT/$PREFIX.1.clean.unpaired.fastq.gz
        CLEAN_UNPAIRED_READ2=$OUTPUT/$PREFIX.2.clean.unpaired.fastq.gz
        THREADS=$2
        # ensure output directory exists
        mkdir -p $OUTPUT

        # trimmomatic command
        java -jar $HOME/jar/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $THREADS $READ1 $READ2 $CLEAN_PAIRED_READ1 $CLEAN_UNPAIRED_READ1 $CLEAN_PAIRED_READ2 $CLEAN_UNPAIRED_READ2 \
        ILLUMINACLIP:$HOME/jar/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:60

#. We need to make that executable (``chmod +x trimmomatic-sub.sh``)
#. Finally, we need to build the submission script to call GNU Parallel and run our job on multiple nodes.  You will need to decide how many nodes to use to trim your samples.  This should parallelize across all nodes.  Be sure to edit the ``CORES_PER_JOB`` to be equivalent to the value you selection above.  Here, we're using 4 cores per job across 2 nodes of 20 cores each (so 5 jobs per node are running simultaneously):

    .. code-block:: bash

        #!/bin/bash       
        #PBS -A <allocation>
        #PBS -l nodes=2:ppn=20
        #PBS -l walltime=2:00:00
        #PBS -q checkpt
        #PBS -N multi_trimmomatic

        # SET THE NUMBER of Cores per job (needs to be multiple of 2)
        export CORES_PER_JOB=4

        # DONT EDIT BELOW #

        # We need java to run trimmomatic
        module load jdk/1.8.0_161
        module load gnuparallel/20170122

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # automatically set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --colsep '\,' \
                --progress \
                --joblog logfile.trimmomatic.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a files-to-trim.txt \
                ./trimmomatic-sub.sh {$1} {$2}
    
#. Once we have trimmed our reads, it's time to setup the genome to which we want to map our reads.  Create a folder (``reference``) to hold the genome, then place the FASTA file for the assembly in this folder.  So, our working directory looks like this:

    .. code-block:: bash

        .
        ├── clean-reads
        ├── files-to-trim.txt
        ├── logfile.629193.smic3
        ├── raw-reads
        ├── reference
        │   └── Tcae_B94639.pseudo.upper.masked.fasta
        ├── trimmomatic.qsub
        └── trimmomatic-sub.sh

#. Now, we need to generate the ``bwa`` index of the genome, as well as the sequence dictionary and fasta index that we'll need for GATK later.  Create ``bwa-index.qsub`` and then submit to the cluster:

    .. code-block:: bash

        #PBS -A <allocation>
        #PBS -l nodes=1:ppn=1
        #PBS -l walltime=12:00:00
        #PBS -q single
        #PBS -N bwa_index

        # ADD THE PATH TO YOUR REFERENCE
        REFERENCE=./reference/Tcae_B94639.pseudo.upper.masked.fasta

        # DONT EDIT BELOW #

        # be sure to activate our conda environment (and we have to use "source"
        # instead of "conda") as well as load java
        source activate mapping
        module load jdk/1.8.0_161

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # create a samtools index
        bwa index $REFERENCE
        samtools faidx $REFERENCE

        # create a sequence dictionary, which we'll need later.  Use full path to GATK because it's 
        # in a different conda environment
        /project/brant/home/miniconda3/envs/gatk/bin/gatk CreateSequenceDictionary --REFERENCE $REFERENCE


#. Once the index is built, we can generate alignment BAMs for each set of reads in the data set.  We'll do that using a similar approach as before.  First, we need to generate a list of files to process along with some associated metadata. Here, it might be easier to use something like excel to generate our sample list - we want 4 columns of information - THREADS, READ1, READ2, REFERENCE, NAME where NAME is the name that you will use to identify each sample in your analysis (analyses).  You can get a head start on the information that you need by running the following.  The last thing you need to do is to add a column of values for the NAME of each sample:
    
    .. code-block:: bash
   
        cd <path to where you are working>
        REFERENCE=$PWD/reference/Tcae_B94639.pseudo.upper.masked.fasta
        ls -d $PWD/clean-reads/*.1.clean.paired.fastq.gz | sed -E "s/(.*).1.clean.paired.fastq.gz/5,\1.1.clean.paired.fastq.gz,\1.2.clean.paired.fastq.gz/" | sed -E "s|.*|&,$REFERENCE|"

#. I have edited the output of the above to add sample names in the last column (it's comma-delimited) and named the file ``files-to-align.txt``.  That file looks like the following:

    .. code-block:: bash

        5,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGCTTT.unmapped.1.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGCTTT.unmapped.2.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/reference/Tcae_B94639.pseudo.upper.masked.fasta,bob
        5,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGGAAT.unmapped.1.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGGAAT.unmapped.2.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/reference/Tcae_B94639.pseudo.upper.masked.fasta,john
        5,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGTGCC.unmapped.1.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGTGCC.unmapped.2.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/reference/Tcae_B94639.pseudo.upper.masked.fasta,steve
        5,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGTTCC.unmapped.1.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.AGTTCC.unmapped.2.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/reference/Tcae_B94639.pseudo.upper.masked.fasta,sue
        5,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.ATCCGC.unmapped.1.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.ATCCGC.unmapped.2.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/reference/Tcae_B94639.pseudo.upper.masked.fasta,sally
        5,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.ATGACT.unmapped.1.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/clean-reads/HC2HMDSXX.1.ATGACT.unmapped.2.clean.paired.fastq.gz,/home/brant/work/tmp/rafa/reference/Tcae_B94639.pseudo.upper.masked.fasta,sarah

#. We need to setup a script we'll call w/ ``Parallel`` to run ``bwa`` and ``samtools`` against each sample to generate a BAM while also adding RG header info to each sample as we align.  Note that in the following, we're using a separate conda environment containing ``bwa`` and ``samtools``:

    .. code-block:: bash

        #!/bin/bash

        # name of output folder
        OUTPUT=bwa-alignments

        ## DO NOT EDIT BELOW THIS LINE - this comes as input from GNU parallel on STDIN ##

        source activate mapping

        THREADS=$1
        READ1=$2
        READ2=$3
        REFERENCE=$4
        SAMPLE_NAME=$5

        ## Here are the specific commands we are running

        # ensure that the output directory exists
        mkdir -p $OUTPUT && cd $OUTPUT

        # create the RG header for each sample to it gets in there while mapping
        # this assume all data are from the same read-group (library)
        HEADER=`printf @RG%sID:%s%sSM:%s%sPL:ILLUMINA '\\t' $SAMPLE_NAME '\\t' $SAMPLE_NAME '\\t'`
        
        # run bwa and output BAM, sort that BAM and index it
        bwa mem -t $THREADS -R "${HEADER}" $REFERENCE $READ1 $READ2 | samtools view -bS - > $SAMPLE_NAME.bam &&
        samtools sort -@ $THREADS $SAMPLE_NAME.bam -o $SAMPLE_NAME.sorted.bam &&
        samtools index -@ $THREADS $SAMPLE_NAME.sorted.bam &&
        rm $SAMPLE_NAME.bam

#. Now, setup the qsub that we need to run the jobs in parallel.  Be sure to enter the number of threads you selected for each job at the top of the script

    .. code-block:: bash

        #!/bin/bash       
        #PBS -A <allocation>
        #PBS -l nodes=2:ppn=20
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N multi_bwa

        # SET THE NUMBER of Cores per job (the number of cores on a node needs to be divisible by this #)
        export CORES_PER_JOB=5

        # DONT EDIT BELOW #

        # load GNU parallel
        module load gnuparallel/20170122

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # automatically set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --colsep '\,' \
                --progress \
                --joblog logfile.align.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a files-to-align.txt \
                ./bwa-align-sub.sh {$1} {$2} {$3} {$4} {$5} 

#. Once that's all finished, your directory structure should look something like this:

    .. code-block:: bash

        .
        ├── bwa-alignments
        │   ├── bob.sorted.bam
        │   ├── bob.sorted.bam.bai
        │   ├── john.sorted.bam
        │   ├── john.sorted.bam.bai
        │   ├── sally.sorted.bam
        │   ├── sally.sorted.bam.bai
        │   ├── sarah.sorted.bam
        │   ├── sarah.sorted.bam.bai
        │   ├── steve.sorted.bam
        │   ├── steve.sorted.bam.bai
        │   ├── sue.sorted.bam
        │   └── sue.sorted.bam.bai
        ├── bwa-align-sub.sh
        ├── bwa_index.e629199
        ├── bwa_index.e629202
        ├── bwa_index.o629199
        ├── bwa_index.o629202
        ├── bwa-index.sh
        ├── bwa.qsub
        ├── clean-reads
        ├── files-to-align.txt
        ├── files-to-trim.txt
        ├── logfile.629193.smic3
        ├── logfile.align.629214.smic3
        ├── multi_trimmomatic.e629193
        ├── multi_trimmomatic.e629208
        ├── multi_trimmomatic.o629193
        ├── multi_trimmomatic.o629208
        ├── raw-reads
        ├── reference
        ├── trimmomatic.qsub
        └── trimmomatic-sub.sh
    
#. Now, we're read to mark duplicates in all of the BAM files. Again, we'll take the same approach as above, although this time we can use a shorter format. Let's create a file of threads, reference, and input file names:

    .. code-block:: bash
    
        cd <path to where you are working>
        REFERENCE=$PWD/reference/Tcae_B94639.pseudo.upper.masked.fasta
        THREADS=5
        for BAM in bwa-alignments/*.sorted.bam; do
            echo $THREADS,$REFERENCE,$PWD/$BAM >> bams-to-clean.txt;
        done


#. Create the script that were going to call with ``Parallel``, save it as ``mark-and-fix-dupes-sub.sh``, and make it executable with ``chmod +x mark-and-fix-dupes-sub.sh``:

    .. code-block:: bash

        #!/bin/bash

        # name of output folder
        OUTPUT=md-alignments

        ## DO NOT EDIT BELOW THIS LINE - this comes as input from GNU parallel on STDIN ##
        module load jdk/1.8.0_161
        source activate gatk

        THREADS=$1
        REFERENCE=$2
        INPUT=$3
        FILENAME=$(basename -- "$INPUT")
        FILENAME_PART="${FILENAME%.*}"
        OUT1=$OUTPUT/$FILENAME_PART.md.bam
        OUT2=$OUTPUT/$FILENAME_PART.md.fx.bam

        # ensure that the output directory exists
        mkdir -p $OUTPUT 
        # run duplicate marking using Spark (in local mode) and setting the cores appropriately
        # we are assuming the number of threads here will be 5
        gatk --java-options "-Xmx16G" MarkDuplicatesSpark --spark-runner LOCAL --input $INPUT --output $OUT1 --conf 'spark.executor.cores=$THREADS' &&
        gatk --java-options "-Xmx16G" SetNmMdAndUqTags --INPUT $OUT1 --OUTPUT $OUT2 --REFERENCE_SEQUENCE $REFERENCE &&
        rm $OUT1 && rm $OUT1.bai && rm $OUT1.sbi &&
        gatk --java-options "-Xmx16G" BuildBamIndex --INPUT $OUT2


#. Setup the qsub script to submit this with GNU Parallel:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <allocation>
        #PBS -l nodes=2:ppn=20
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N multi_mark_dupe

        # SET THE NUMBER of Cores per job (the number of cores on a node needs to be divisible by this #)
        export CORES_PER_JOB=5

        # DONT EDIT BELOW #

        # load GNU parallel
        module load gnuparallel/20170122

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # automatically set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --colsep '\,' \
                --progress \
                --joblog logfile.dupes.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a bams-to-clean.txt \
                ./mark-and-fix-dupes-sub.sh {$1} {$2} {$3}

#. If you already have set of very high-quality SNPs that you can perform Variant Quality Score Recalibration (VQSR) with at this stage - do that (see below).  We will assume that you do not have these, so you need to go through at least one round of SNP calling to generate this set. That begins with ``HaplotypeCaller`` in GVCF-output mode, which we will run in single-threads, but setting the RAM for each thread file to 4 GB.  On @smic, this means we can run a total of 16 threads.  First, make file that contains the path to our REFERENCE sequence and each BAM file:

    .. code-block:: bash

        REFERENCE=$PWD/reference/Tcae_B94639.pseudo.upper.masked.fasta
        for BAM in md-alignments/*.md.fx.bam; do 
            echo "$REFERENCE,$BAM" >> bams-to-haplotype-call.txt;
        done

#. Now, setup the script that GNU Parallel will call, save it as ``haplotype-gvcf-sub.sh``, and make it executable ``chmod +x haplotype-gvcf-sub.sh``:

    .. code-block:: bash

        #!/bin/bash

        # name of output folder
        OUTPUT=temp-gvcf

        ## DO NOT EDIT BELOW THIS LINE - this comes as input from GNU parallel on STDIN ##
        module load jdk/1.8.0_161
        source activate gatk

        REFERENCE=$1
        INPUT=$2
        FILENAME=$(basename -- "$INPUT")
        FILENAME_PART="${FILENAME%.*}"
        OUT1=$OUTPUT/$FILENAME.g.vcf.gz 

        # ensure that the output directory exists
        mkdir -p $OUTPUT 
        # NOTE - we're assuming single-threaded operation for each BAM file, so we set RAM to 4GB each (16 cores, max)
        gatk --java-options "-Xmx4G" HaplotypeCaller -R $REFERENCE -I $INPUT -O $OUT1 -ERC GVCF

#. Finally, setup the GNU parallel qsub script.  Each node in the following will run 16 BAMs:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A <allocation>
        #PBS -l nodes=1:ppn=20
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N multi_haplotype_call

        # SET THE NUMBER of Cores per job (the number of cores on a node needs to be divisible by this #)
        export CORES_PER_JOB=1

        # DONT EDIT BELOW #

        # load GNU parallel
        module load gnuparallel/20170122

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # automatically set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --colsep '\,' \
                --progress \
                --joblog logfile.haplotype_gvcf.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a bams-to-haplotype-call.txt \
                ./haplotype-gvcf-sub.sh {$1} {$2}

#. Now, we need to integrate the GVCF files together, and the first step of that uses ``GenomicsDBImport``. ``GenomicsDBImport`` requires a tab delimited list of sample names and file names for each sample, so generate that:

    .. code-block:: bash

        for VCF in temp-gvcf/*; do
            filename=$(basename -- "$VCF");
            name="${filename%%.*}";
            echo -e "$name\t$VCF" >> gvcfs-for-db-import.sample_map;
        done

#. Now we can run ``GenomicsDBImport`` to bring together all the gvcf files.  You will run this normally using @smic.  This is not multithreaded:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A hpc_allbirds04
        #PBS -l nodes=1:ppn=20
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N GenomicsDBImport

        # activate gatk
        source activate gatk

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # we need a tmp dir that is large for the program to use
        mkdir -p $PBS_O_WORKDIR/tmp

        # set batch size equal to cores on node also set max RAM a 
        # little low, because there is additional overhead
        # involved per GATK website
        gatk --java-options "-Xmx58g" \
            GenomicsDBImport \
            --genomicsdb-workspace-path my_database \
            --tmp-dir=$PBS_O_WORKDIR/tmp \
            --batch-size 20 \
            --sample-name-map gvcfs-for-db-import.sample_map

#. Once that is run, we will call the population of genotypes in all samples using ``GenotypeGVCFs`` run against the database of all individuals:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A hpc_allbirds04
        #PBS -l nodes=1:ppn=20
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N GenotypeGVCFs

        # activate gatk
        source activate gatk

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # set reference
        REFERENCE=$PWD/reference/Tcae_B94639.pseudo.upper.masked.fasta

        # now go ahead and run and set batch size equal to cores on node
        # also set max RAM a little low, because there is additional overhead
        # involved
        gatk --java-options "-Xmx58g" \
            GenotypeGVCFs \
            -R $REFERENCE \
            -V gendb://my_database \
            -O output.vcf.gz

#. This will output a file of genotypes for all individuals that you will need to statically filter to keep only the best genotypes.  You can do this using ``vcftools``, and we have recently used a command similar to the following to identify the "best SNPs" for BQSR.  Depending on your data set, you might want to adjust some of these parameters (in particular, ``--max-missing``).  You might also consider additional filtering options based on your particular data set (e.g. are loci in HWE).

    .. code-block:: bash

        vcftools \
            --vcf output.vcf.gz \
            --minDP 30 \
            --minQ 30 \
            --minGQ 30 \
            --max-alleles 2 \
            --remove-indels \
            --max-missing 0.5

#. Now, you should review the GATK article on BQSR.  We can use the valid SNPs to perform BQSR, and we need to return to our original BAM files because these are what we are recalibrating.  First thing we need to do is to make an input file listing the REFERENCE, the ``--known-sites``, and the BAM:

    .. code-block:: bash

        REFERENCE=$PWD/reference/Tcae_B94639.pseudo.upper.masked.fasta
        SITES=$PWD/best.output.vcf.gz
        for BAM in md-alignments/*; do
            echo "$REFERENCE,$SITES,$BAM" >> bams-to-recalibrate.txt;
        done

#. Now, create a script we'll run w/ GNU parallel which implements the first and second stage of BQSR for each BAM file.  Again, we'll limit the RAM for each BAM to 4 GB, so we can run 16 files in parallel.  We will name this file ``bam-bqsr-sub.sh``, and we need to ``chmod +x bam-bqsr-sub.sh`` after creating the file with the following contents:

    .. code-block:: bash

        #!/bin/bash

        # name of output folder
        OUTPUT=bqsr-alignments

        ## DO NOT EDIT BELOW THIS LINE - this comes as input from GNU parallel on STDIN ##
        module load jdk/1.8.0_161
        source activate gatk

        REFERENCE=$1
        SITES=$2
        INPUT=$3
        FILENAME=$(basename -- "$INPUT")
        FILENAME_PART="${FILENAME%.*}"
        OUT1=$OUTPUT/$FILENAME.recal_data.table
        OUT2=$OUTPUT/$FILENAME.bqsr.bam

        # ensure that the output directory exists
        mkdir -p $OUTPUT 
        # NOTE - we're assuming single-threaded operation for each BAM file, 
        # so we set RAM to 4GB each (16 cores, max)
        gatk --java-options "-Xmx4G" BaseRecalibrator \
            -I $INPUT \
            -R $REFERENCE \
            --known-sites $SITES \
            -O $OUT1 \
        && gatk --java-options "-Xmx4G" ApplyBQSR \
            -R $REFERENCE \
            -I $INPUT \
            --bqsr-recal-file $OUT1 \
            -O $OUT2

#. Setup the GNU Parallel script to run the ``bam-bqsr-sub.sh``:

    .. code-block:: bash

        #!/bin/bash
        #PBS -A hpc_allbirds04
        #PBS -l nodes=1:ppn=20
        #PBS -l walltime=4:00:00
        #PBS -q checkpt
        #PBS -N multi_bqsr

        # SET THE NUMBER of Cores per job (the number of cores on a node needs to be divisible by this #)
        export CORES_PER_JOB=1

        # DONT EDIT BELOW #

        # load GNU parallel
        module load gnuparallel/20170122

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # automatically set the number of Jobs per node based on $CORES_PER_JOB
        export JOBS_PER_NODE=$(($PBS_NUM_PPN / $CORES_PER_JOB))

        parallel --colsep '\,' \
                --progress \
                --joblog logfile.bqsr.$PBS_JOBID \
                -j $JOBS_PER_NODE \
                --slf $PBS_NODEFILE \
                --workdir $PBS_O_WORKDIR \
                -a bams-to-recalibrate.txt \
                ./bam-bqsr-sub.sh {$1} {$2} {$3}