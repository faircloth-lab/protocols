.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Calling SNPs with GATK
======================

:Author: Jessie Salter, Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Calling SNPs with GATK`_ 

.. _Calling SNPs with GATK: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-gatk.rst

Purpose
-------

Often, for phylogeneographic and population genetics studies, we want to call SNPs.  There are several ways to do that, but the following protocol is meant to document how we call and filter SNPs with GATK_.

Preliminary Steps
-----------------

#. You will want to "clean" your raw reads by trimming them for adapter contamination and low quality reads.  This is pretty standard and can easily be accomplished for lots of reads using trimmomatic_.

#. You will also want a genome/reference sequence to which you will be aligning raw reads.  It's best to have repeats in this assembly soft-masked with something like repeatmasker_.  Most genomes that come from sources like UCSC and NCBI area already repeat masked. New assemblies are not, and generally should be repeat masked.

#. I setup a conda environment just for SNP calling.  I did that with the following:

    .. code-block:: bash
    
        # create an environment
        conda create -n gatk python=3.7 bwa samtools vcftools openjdk=8.0.192
        # activate
        conda activate gatk
        # gatk has annoying distribution rights, so we need to get our
        # own copy and unpack that.  we'll keep it in /<conda-env>/opt
        mkdir /project/brant/home/miniconda3/envs/gatk/opt
        cd /project/brant/home/miniconda3/envs/gatk/opt
        wget https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip
        unzip gatk-4.1.7.0.zip
        # you can symlink into /<conda-env>/bin for ease of use
        cd ../bin
        ln -s ../opt/gatk-4.1.7.0/gatk
        # once environment is activated, you should be able to run `gatk`
        gatk --version
        
Steps
-----

#. The first step of the process is to get your raw read data and trim those for adapter contamination and low quality bases.  Although this is not required as part of GATK, our read data are sometimes short and have adapter read-through at the 3' end.  So, it's best to trim, and it's probably easiest to trim your data using _illumiprocessor.  That's already been done in the example, below.

#. For the example, we are aligning data to a single chromosome of an avian assembly.  I've generally done that to make things reasonable fast/easy in terms of testing.  Those data are on @smic.  The reference contig is in ``reference`` and the raw read data are in ``trimmed-fastq``:

    .. code-block:: bash

        .
        ├── reference
        │   └── ScNXozI_1979.fasta
        └── trimmed-fastq
            ├── LSU_11071
            │   ├── LSU_11071.ScNXozI_1979.r1.fastq.gz
            │   └── LSU_11071.ScNXozI_1979.r2.fastq.gz
            ├── LSU_73043
            │   ├── LSU_73043.ScNXozI_1979.r1.fastq.gz
            │   └── LSU_73043.ScNXozI_1979.r2.fastq.gz
            ├── LSU_80496
            │   ├── LSU_80496.ScNXozI_1979.r1.fastq.gz
            │   └── LSU_80496.ScNXozI_1979.r2.fastq.gz
            ├── LSU_80711
            │   ├── LSU_80711.ScNXozI_1979.r1.fastq.gz
            │   └── LSU_80711.ScNXozI_1979.r2.fastq.gz
            └── LSU_81131
                ├── LSU_81131.ScNXozI_1979.r1.fastq.gz
                └── LSU_81131.ScNXozI_1979.r2.fastq.gz

#. Before moving on too much, we need to prepare our reference "genome" for alignment.  This usually requires a few steps with _BWA and _Picard - we'll go ahead and run both.  We'll put this script in ``reference``, and it will run from that directory.  We'll run this in the ``single`` queue because it's not multithreaded (and shouldn't gobble much RAM). With larger genomes, this can take a little while to run, but once you've indexed a genome, you don't generally need to do it again, as long as the assembly stays the same.

    .. code-block:: bash

        #PBS -A hpc_allbirds03
        #PBS -l nodes=1:ppn=1
        #PBS -l walltime=12:00:00
        #PBS -q single
        #PBS -N bwa_index

        # be sure to activate our conda environment (and we have to use "source"
        # instead of "conda")
        source activate gatk

        # move into the directory containing this script
        cd $PBS_O_WORKDIR

        # create a samtools index
        bwa index ScNXozI_1979.fasta

        # create a sequence dictionary, which we'll need later
        gatk CreateSequenceDictionary --REFERENCE ScNXozI_1979.fasta

#. Once that's finished, we need to map our raw read data to the reference assembly that we just indexed.  How you pick the reference to map to is a matter of some debate (and not covered here). Keep in mind that once you make a selection, it's rather hard to move between reference assemblies (at least, it's always harder than you would like for it to be).  We usually use _BWA as our mapping program of choice.  After indexing our assembly and trimming our read, our data folders look something like this:

    .. code-block:: bash

.
├── reference
│   ├── bwa_index.e562696
│   ├── bwa_index.o562696
│   ├── bwa_index.qsub
│   ├── ScNXozI_1979.dict
│   ├── ScNXozI_1979.fasta
│   ├── ScNXozI_1979.fasta.amb
│   ├── ScNXozI_1979.fasta.ann
│   ├── ScNXozI_1979.fasta.bwt
│   ├── ScNXozI_1979.fasta.pac
│   └── ScNXozI_1979.fasta.sa
└── trimmed-fastq
    ├── LSU_11071
    │   ├── LSU_11071.ScNXozI_1979.r1.fastq.gz
    │   └── LSU_11071.ScNXozI_1979.r2.fastq.gz
    ├── LSU_73043
    │   ├── LSU_73043.ScNXozI_1979.r1.fastq.gz
    │   └── LSU_73043.ScNXozI_1979.r2.fastq.gz
    ├── LSU_80496
    │   ├── LSU_80496.ScNXozI_1979.r1.fastq.gz
    │   └── LSU_80496.ScNXozI_1979.r2.fastq.gz
    ├── LSU_80711
    │   ├── LSU_80711.ScNXozI_1979.r1.fastq.gz
    │   └── LSU_80711.ScNXozI_1979.r2.fastq.gz
    └── LSU_81131
        ├── LSU_81131.ScNXozI_1979.r1.fastq.gz
        └── LSU_81131.ScNXozI_1979.r2.fastq.gz

If your data are in a format similar to those output by trimmomatic (reads within names directories - like the above), you can map across those directories of raw reads in serial, using something like the following (see below for parallel execution).  This scripts assumes you are using a single node with 20 compute cores, and that your directory structure is similar to the above.  

    .. admonition:: Note
    
    Note that in the following, we're adding the ``RG`` ("read group) headers to each BAM file we create - this is basically a way that downstream programs can identify "individual-specific" data when we do things like combine BAM files.  For more information on RG headers, see `this document <https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups>`_.

#. Before running, create a directory named ``bams``, and this is where we'll be working from.  Copy the script below into ``bams`` as ``bwa_align and submit it with qsub.

.. code-block:: bash

#PBS -A hpc_allbirds03
#PBS -l nodes=1:ppn=20
#PBS -l walltime=12:00:00
#PBS -q checkpt
#PBS -N bwa_align

source activate gatk

# move into the directory containing this script
cd $PBS_O_WORKDIR

for i in ../trimmed-fastq/*;
do
    SAMPL=$(basename $i);
    GENOME=../reference/ScNXozI_1979.fasta;
    READ1=$i/${SAMPL}.ScNXozI_1979.r1.fastq.gz;
    READ2=$i/${SAMPL}.ScNXozI_1979.r1.fastq.gz;
    # note that this is adding the RG header information directly to the SAM/BAM on creation
    bwa mem -t  20 -R `@RG\tID:${SAMPL}\tSM:${SAMPL}` $GENOME $READ1 $READ2 | samtools view -bS - > $SAMPL.bam;
done


for i in ../trimmed-fastq/*;
do
    SAMPL=$(basename $i);
    GENOME=../reference/ScNXozI_1979.fasta;
    READ1=$i/${SAMPL}.ScNXozI_1979.r1.fastq.gz;
    READ2=$i/${SAMPL}.ScNXozI_1979.r1.fastq.gz;
    echo `@RG\tID:${SAMPL}\tSM:${SAMPL}` $GENOME $READ1 $READ2
done


And that our clean reads are in `/path/to/clean/reads/`:






First, get a list of the directories that contain the read data you'll be mapping:


Then, setup your script that will run 

    .. code-block:: bash

        

