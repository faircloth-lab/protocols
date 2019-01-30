.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


Running PAUP
============

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Running PAUP`_ 

.. _Running PAUP: http://github.com/faircloth-lab/protocols/commits/master/protocols-computer/analysis/analysis-paup.rst

Purpose
-------

Often, prior to running more computationally intensive analyses of large phylogenies, we'll take a look at parsimony trees so that we can make sure things are reasonably sensible before moving ahead.  These are pretty simple instructions for generating a parsimony tree using PAUP.


Preliminary Steps
-----------------

#. Before using PAUP, you need to have the PAUP binary availble on your computer.  You can download the binary `from here <http://phylosolutions.com/paup-test/>`_.  You can usually do this using a tool like ``wget``.  If you are running on our local machines, you will want the CentOS X86_64 version.

#. Once you have the binary on your computer, you need to make sure it is in your ``$PATH``, and that it is set to be executable.  Usually, if you place the binary in ``$HOME/bin``, that will be in your path.  Then, you need to ``chmod 0755 <binary name>``.

Steps
-----

#. PAUP requires an input file in NEXUS format.  You can produce this in phyluce using:

    .. code-block:: bash

        phyluce_align_format_nexus_files_for_raxml \
                --alignments mafft-fasta-trimal-clean-75p-complete \
                --output mafft-fasta-trimal-clean-75p-complete-nexus \
                --nexus

#. Once you have a NEXUS-formatted input file, you need to start PAUP on the command line.  Assuming PAUP is in your ``$PATH``, is executable, and is named ``paup``, run:

    .. code-block:: bash

        paup

#. Then, read in the NEXUS formatted alignment file using:

    .. code-block:: bash

        execute name-of-your-file.nexus;

#. Set the parsimony criterion, tell PAUP to root on the outgroup, set the outgroup taxon (here, replacing  <name_of_tip> with an actualy tip in your tree/alignment):

    .. code-block:: bash

        set criterion=parsimony;
        set root=outgroup;
        set storebrlens=yes;
        set increase=auto;
        outgroup <name_of_tip>;

#. Now, run the search; save the results to a file, replacing ``<output_file_name>`` with the output tree name you want; and quit:

    .. code-block:: bash

        hsearch mulTrees=No;
        savetrees file=<output_file_name>.tre format=altnex brlens=yes;
        quit;



