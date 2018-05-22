.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _BioSamples:

Register NCBI BioSamples for a BioProject
=========================================

:Author: Carl Oliveros, Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Register NCBI BioSamples for a BioProject History`_ 

.. _Register NCBI BioSamples for a BioProject History: http://github.com/faircloth-lab/protocols/commits/master/protocols-data/ncbi/register-a-biosample.rst


Purpose
-------

Prior to submitting massively parallel sequencing (MPS) data to NCBI, you want to register (1) your project and (2) those samples involved in a given project.  Registering your samples as [NCBI BioSamples](https://www.ncbi.nlm.nih.gov/biosample/) is the second step of this process (the first is to [[Register an NCBI BioProject]]).  You can register a BioProject and BioSamples before you plan to upload your data.  That said, our laboratory policy is to **upload all the data, all the time** prior to publication (and usually prior to making a paper available as a pre-print).

Steps
-----

#. If you have not done so, :ref:`BioProject`.

#. Log in to the `NCBI submissions portal <https://submit.ncbi.nlm.nih.gov/>`_, and click on the BioSample link:

    .. image:: /images/ncbi-biosample.png

#. Click on the `New submission` link. The **Submitter** info should be the same as what you used for your BioProject.  Click `Continue`.

#. On the **General Info** tab, specify a `Release Date` and select `Batch/Multiple BioSamples`  Click `Continue`.

    .. image:: /images/ncbi-biosample-general-info.png

#. On the **Sample type** tab, select `Model organism or animal sample`.  Click `Continue`.

    .. image:: /images/ncbi-biosample-sample-type.png

#. On the **Attributes** tab, download the Excel template for your BioSample.

    .. image:: /images/ncbi-biosample-attributes.png

    Follow the instructions on the worksheet.  Most of the fields have a detailed description if you hover your pointer over the top right corner of the cell.

    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Field                                                     | Comments                                                                                                                                                                    |
    +===========================================================+=============================================================================================================================================================================+
    | sample_name                                               | Use something like "Species-genus-Institution-Accession"                                                                                                                    |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | BioProject accession                                      | There should be only one for all rows                                                                                                                                       |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | organism                                                  | e.g. "Gallus gallus", as much as possible, a valid NCBI taxonomy name                                                                                                       |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | strain, isolate, breed, cultivar, ecotype, age, dev_stage | You can fill these with "not applicable" or "not collected"                                                                                                                 |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | sex                                                       | “male", "female", or "not collected"                                                                                                                                        |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | tissue                                                    | muscle, liver, toe pad                                                                                                                                                      |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | geo_loc_name                                              | tissue or toe pad                                                                                                                                                           |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | specimen_voucher                                          | Use "Institution:Accession", e.g. "LSUMZ:Ornithology:12345". A list of valid institution codes can be found `here <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/coll_dump.txt>`_ |
    +-----------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

    .. note:: You should have one BioSample for each specimen, and each of your BioSamples must have differentiating information (excluding sample name, title, bioproject accession and description). This check was implemented to encourage submitters to include distinguishing information in their samples. If the distinguishing information is in the sample name, title or description, please recode it into an appropriate attribute, either one of the predefined attributes or a custom attribute you define. If it is necessary to represent true biological replicates as separate BioSamples, you might add an 'aliquot' or 'replicate' attribute, e.g., 'replicate = biological replicate 1', as appropriate.

#. Once the Excel sheet is completed, save it as a tab delimited text file and upload it on the **Attributes** tab.  Click `Continue`.

   The **Comments** tab will tell you any errors or warnings associated with your BioSample worksheet.  Be sure to correct all errors before continuing.  Warnings may be ok.  Click `Continue`.  Review the submission on the Overview tab and click `Submit` at the bottom if no other changes are necessary.  Wait for the BioSample to be processed.

   Sometimes, taxonomy differences between your "organism" (from table above) can conflict with the entries in `NCBI Taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_, and sometimes the changes needed can get hung up at NCBI.  If you've been waiting for more than 3-4 business days for your BioSamples to process, you should email NCBI at the contact for BioSamples (biosamplehelp@ncbi.nlm.nih.gov).


#. You should eventually receive an email from NCBI that looks similar to this:
 
    .. code-block:: text

        Dear Brant Faircloth,

        This is an automatic acknowledgment that your recent submission to the BioSample
        database has been successfully processed and will be released on the date
        specified.

        BioSample accessions:  SAMN05915021, SAMN05915022, ... see attached file.
        Temporary SubmissionID:  SUB2020739
        Release date:  when referenced data is published

        Your BioSample records will be accessible with the following links:
        See object links in the attachment to this message.

        Please reference BioSample accessions SAMN05915021, SAMN05915022, SAMN05915023,
        SAMN05915024, SAMN05915025, SAMN05915026, SAMN05915027, SAMN05915028,
        SAMN05915029, SAMN05915030, ... see attached file. when making corresponding
        sequence data submissions.

        Send questions and update requests to biosamplehelp@ncbi.nlm.nih.gov; include
        the BioSample accessions SAMN05915021, SAMN05915022, SAMN05915023, 
        SAMN05915024, SAMN05915025, SAMN05915026, SAMN05915027, SAMN05915028,
        SAMN05915029, SAMN05915030, ... see attached file. in any correspondence.

        Regards,

        NCBI BioSample Submissions Staff
        Bethesda, Maryland USA
        ***********************************************************
        (301) 496-2475
        (301) 480-2918 (Fax)
        biosamplehelp@ncbi.nlm.nih.gov (for BioSample questions/replies)
        info@ncbi.nlm.nih.gov (for general questions regarding NCBI)
        ***********************************************************


#. If you are uploading data to NCBI SRA proceed to :ref:`Submitting Read Data to NCBI SRA`_.