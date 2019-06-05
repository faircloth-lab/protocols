.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _ArchiveSequencingData:

Archiving Sequencing Data
=========================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Archiving Sequencing Data History`_.

.. _Archiving Sequencing Data History: https://github.com/faircloth-lab/protocols/commits/master/protocols-data/lab/archive-sequencing-data.rst


Purpose
-------

We archive the raw read data from ALL of our old sequencing runs (>5 years) in AWS Glacier Deep Archive.


Steps
-----

#. Navigate to the ``illumina-runs`` section of the NFS
#. Identify sequencing directories needed to be archived.  Usually, in ``illumina-runs``, this is everything but the ``clean`` sequence directory that is ~5 years old. Once you've decided what needs to be excluded for a particular sequencing run, generate a directory tree of what you're uploading:

    .. code-block:: bash

        tree -I 'clean' -a > ${PWD##*/}-dirtree.txt

#. Check to make sure all files in dirs to package up are already zipped (these are sequence files, so they should be)
#. Now, package up everything, except for any excluded directories (e.g. like ``clean``, which just duplicates the data):

    .. code-block:: bash

        tar --exclude ./clean -cvf ${PWD##*/}.tar ./ | tee ${PWD##*/}-tar.out

#. This will make an output file named ``<directory-name>-tar.out`` that you can check to ensure everything has been packaged up that you wanted
#. Compute md5 checksums of everything:

    .. code-block:: bash

        md5sum ${PWD##*/}.tar ${PWD##*/}-tar.out ${PWD##*/}-dirtree.txt > ${PWD##*/}.md5

#. Now, go ahead and upload those to AWS Glacier Deep Archive (be sure to use the correct ``--profile``:

    .. code-block:: bash
    
        aws s3 cp ${PWD##*/}.tar s3://2013-faircloth-lab-sequence-data --storage-class DEEP_ARCHIVE --profile lab-data &&  
        aws s3 cp ${PWD##*/}-tar.out s3://2013-faircloth-lab-sequence-data --storage-class DEEP_ARCHIVE --profile lab-data &&
        aws s3 cp ${PWD##*/}-dirtree.txt s3://2013-faircloth-lab-sequence-data --storage-class DEEP_ARCHIVE --profile lab-data &&
        aws s3 cp ${PWD##*/}.md5 s3://2013-faircloth-lab-sequence-data --storage-class DEEP_ARCHIVE --profile lab-data
        
#. Check the tar log file to make sure everything got packaged up (basically avoiding the ``clean`` data).
#. Remove the directories that you've archived.  For now, leave the ``clean`` data directory
#. Denote in the Google Sheet that the data have been archived
#. Move the directory to ``/nfs/data1/illumina-runs/ARCHIVED``  