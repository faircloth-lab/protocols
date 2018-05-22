.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _BioProject:

Register an NCBI BioProject
===========================

:Author: Brant C. Faircloth, Carl Oliveros
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Register an NCBI BioProject History`_.

.. _Register an NCBI BioProject History: https://github.com/faircloth-lab/protocols/commits/master/protocols-data/ncbi/register-a-bioproject.rst


Purpose
-------

Prior to submitting massively parallel sequencing (MPS) data to NCBI, you want to register (1) your project and (2) those samples involved in a given project.  Registering an `NCBI BioProject`_. is the first step of this process.  You can also register a BioProject before you plan to register your samples or upload your data.  That said, our laboratory policy is to **upload all the data, all the time** prior to publication (and usually prior to making a paper available as a pre-print).

Steps
-----

#. Create an NCBI account at the following link, if you do not have one, at the `NCBI website <https://www.ncbi.nlm.nih.gov/account/register/>`_

#. Log in to the `NCBI submissions portal <https://submit.ncbi.nlm.nih.gov/>`_, and click on the BioProject link
   
    .. image:: /images/ncbi-bioproject.png
    
#. Click on the `New submission` link. The `Submitter` info should automatically be filled in using info from your profile.  Verify the info and click `Continue`.

    .. image:: /images/ncbi-bioproject-submitter.png

#. On the **Project type** tab, select `Raw sequence reads` under **Project data type** and `Multispecies` under **Sample scope**.  Click `Continue`.

    .. image:: /images/ncbi-bioproject-type.png

#. On the **Target** tab, fill in descriptions for `Organism name` and `Multispecies description`.  Click `Continue`.  `Organism name` can be a general description of the class of organisms (e.g. "Oscine passerines").

    .. image:: /images/ncbi-bioproject-target-tab.png

#. On the **General Info** tab, provide a `Release Date`, a `Project Title` (which could be your manuscript title), and `Public Description` (which could be your manuscript abstract).  Select `Evolution` under **Relevance**.  Under this tab, you can also add grant numbers and titles for grants that contributed to the project. If you are supported by NSF or NIH, please add your grant numbers and titles here. Click `Continue`.

     .. image:: /images/ncbi-bioproject-general-info.png

#. On the **BioSample** tab, do not enter any information (we will deal with this for multiple samples as part of [[Register NCBI BioSamples for a BioProject]]) and click `Continue`.

    .. image:: /images/ncbi-bioproject-biosample.png

#. On the **Publications** tab, provide the DOI of your manuscript, if available.  Click `Continue`.

    .. image:: /images/ncbi-bioproject-publications.png

#. Review the information on the **Overview** tab and click `Submit` at the bottom, if no changes are needed.  Wait for the BioProject submission to be processed (may take only a few minutes to receive the email from NCBI).  A processed BioProject will look like this:

    .. image:: /images/ncbi-bioproject-processed.png

#. And, you'll receive an email from NCBI with something similar to the following contents:

    .. code-block:: text

        Dear xxxxx,

        This is an automatic acknowledgment that your submission:

        SubmissionID:       SUB1211501
        BioProject ID:      PRJNA304409
        Title:          

        has been successfully registered with the BioProject database. After review by
        the database staff, your project information will be accessible with the
        following link, usually within a few days of the release date that you set (or
        the release of linked data, whichever is first):

        http://www.ncbi.nlm.nih.gov/bioproject/304409

        Please use the BioProject ID PRJNA304409 with your correspondence and your data
        submissions.

        Send questions to bioprojecthelp@ncbi.nlm.nih.gov, and include the BioProject 
        ID and organism name.

        Regards,

        NCBI BioProject Submissions Staff
        Bethesda, Maryland  USA
        ***********************************************************
        (301) 496-2475
        (301) 480-2918 (Fax)
        bioprojecthelp@ncbi.nlm.nih.gov (for BioProject questions/replies)
        info@ncbi.nlm.nih.gov (for general questions regarding NCBI)
        ***********************************************************

#. If you are registering BioSamples at this time, proceed to :ref:`BioSamples`.  Otherwise, you can use this BioProject number (PRJNAXXXXXXXX) in your manuscript as the "pointer" to all data associated with your project.  A finished BioProject page that points to all available data looks something like this:

    .. image:: /images/ncbi-bioproject-example.png

    This BioProject page provides links to **ALL** NCBI resources related to your BioProject - making the BioProject a one-stop-shop for all of your project sequence data.