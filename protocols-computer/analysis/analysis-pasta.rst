.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running Pasta
=============

:Author: Carl Oliveros and Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running Pasta`_ 

.. _Running Pasta: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-pasta.rst

Purpose
-------

Pasta_ is essentially an updated version of Saté_, and Pasta_ may or may not be incorporated with Treeshrink_ (as we do here).  This combination of programs aligns sequences, removes really bad sequences from those alignments, then puts the cleaned up alignments back through a second alignment round.

.. _Pasta: https://github.com/smirarab/pasta
.. _Saté: http://phylo.bio.ku.edu/software/sate/sate.html
.. _Treeshrink: https://github.com/uym2/TreeShrink

Preliminary Steps
-----------------

#. The following assumes that you've installed phyluce_ correctly (along with ``conda`` and ``bionconda``
#. If you have not already, create a ``conda`` environment for Pasta_ and install Pasta_

    .. code-block:: bash

        conda create -n pasta python=3 pasta

#. Switch to that environment and install treeshrink:

    .. code-block:: bash

        source activate pasta
        conda install -c smirarab treeshrink

#. To leave the environment, run:

    .. code-block:: bash

        source deactivate pasta


Steps
-----

#. Pasta_ (or Pasta_ + TreeShrink_) enter the equation when we're trying to align DNA sequences.  Here, I'll discuss these steps in the context of using the phyluce_ pipeline, although many of the steps in the apparoach are similar regardless of whether you are using phyluce_ or not.

#. You can implement Pasta_ at several stages of the phyluce_ pipeline, but the easiest is probably after you have identified the UCE loci and extracted those loci what we call a "monolithic" fasta file.  First thing you want to do it "explode" that monolithic FASTA:

    .. code-block:: bash

        phyluce_assembly_explode_get_fastas_file \
            --input my-monolithic.fasta \
            --output exploded-loci

#. Now what you should have is a directory of fasta files, one for each locus in your data set.  You likely want to filter these loci to remove really short stuff - typically something like those sequences having < 50% of the median length of all sequences for a particular locus:

    .. code-block:: bash

        python ~/git/phyluce/bin/assembly/phyluce_assembly_filter_seqs_from_fastas \
            --input exploded-loci \
            --output exploded-loci-length-filtered \
            --filtered-sequences-file exploded-loci_fasta.shorts \
            --proportion 0.5 \
            --cores 12 \
            --log-path log

#. Now, we want to filter those FASTA files for a locus that have fewer than 4 taxa.  We can do this using some code meant for alignment data (but adjusting it, here



