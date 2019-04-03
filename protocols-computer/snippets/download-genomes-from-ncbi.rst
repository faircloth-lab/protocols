.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

Download Genomes From NCBI
==========================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

All of the following assume that you are using the Z shell (zsh).  These may or may not work in BASH.


Modification History
--------------------

See `Download Genomes From NCBI`_ 

.. _Download Genomes From NCBI: https://github.com/faircloth-lab/protocols/commits/master/protocols-computer/snippets/download-genomes-from-ncbi.rst

Purpose
-------

Sometimes, we want to download genomes from NCBI and sometimes we need to do that for a lot of genomes.  Thing is, it's not always easy to do this, and NCBI does not make it abundantly clear what the best way to do this is for larger genomes (e.g. not microbes).  So, here's one way to go about it.

Steps
-----

1. You need to find identifier information for the genome(s) you want to download.  NCBI indexes genomes in several ways, some of them weird.  Probably the best way to find the genomes you want is to make a list of the taxa that you want to download (e.g. genus and species).  Then you can feed that list of taxa into the script below to pull down the NCBI Taxonomy ID for that individual taxon.  We'll then use this Taxonomy ID to find the assemblies we're after.  You want your list of taxa to look something like:

    .. code-block:: text

        Benthosema glaciale
        Percopsis transmontana
        Typhlichthys subterraneus
        Cyttopsis rosea
        Gadiculus argenteus
        Trisopterus minutus
        Brosme brosme
        Molva molva
        Phycis phycis
        Phycis blennoides

2. Save that list as 'taxa.txt'. Now, create a new file (``get_tax_id.py``), edit the following to add your email address, and run the following code against this list with Python (the list name is hardcoded into the code below, but it's easy to edit):

    .. code-block:: python

        #!/usr/bin/env python
        # -*- coding: utf-8 -*-

        """
        (c) 2016 Brant Faircloth || http://faircloth-lab.org/
        All rights reserved.

        This code is distributed under a 3-clause BSD license. Please see
        LICENSE.txt for more information.

        Created on 28 April 2016 08:42 CDT (-0500)
        """


        import os
        import sys
        import time
        import argparse

        from Bio import Entrez

        # import pdb

        def get_tax_id(species):
            """to get data from ncbi taxomomy, we need to have the taxid.  we can
            get that by passing the species name to esearch, which will return
            the tax id"""
            species = species.replace(" ", "+").strip()
            search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
            record = Entrez.read(search)
            return record['IdList'][0]


        def get_tax_data(taxid):
            """once we have the taxid, we can fetch the record"""
            search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
            return Entrez.read(search)


        def main():
            Entrez.email = "name@domain.edu"
            if not Entrez.email:
                print("You must add your email address")
                sys.exit(2)
            with open('taxa.txt') as infile:
                all_taxa = {}
                for line in infile:
                    tax_name = line.strip()
                    taxid = get_tax_id(tax_name)
                    print("{},{}".format(tax_name, taxid))
                    time.sleep(1)


        if __name__ == '__main__':
            main()

3. This will spit out a list of taxa to ``stdout`` that looks like:

    .. code-block:: text

        Malacocephalus occidentalis,630739
        Macrourus berglax,473319
        Bathygadus melanobranchus,630650
        Laemonema laureysi,1784819
        Trachyrincus scabrus,562814
        Muraenolepis marmoratus,487677
        Melanonus zugmayeri,181410

4. This list now contains the taxon you searched for, and the NCBI Taxonomy ID for that species.  The code will hit an error if you include a species name that does not exist in the NCBI Taxonomy database.

5. Save the list that's output to a ``csv`` file named something like ``ncbi_id.csv``.  Once you've done that, we need to create a new (potentially temporary) conda_ environment to hold the `NCBI Genome Download <https://github.com/kblin/ncbi-genome-download>`_ code.

    .. code-block:: bash


        conda create -n ncbi python=3 pip
        conda activate ncbi
        pip install ncbi-genome-download

6. With that environment installed, either navigate to (or create) a directory to hold the genomes we want to download, and copy the list output from our automated search against NCBI taxonomy.

7. In this directory, we'll use a ZSH shell script to parse the list of species and NCBI Taxonomy ID we just created, and use components of those parsed files to download the genome sequences we want.  In the example below, we're telling ``ncbi-genome-download`` to download only the assembly-report and the genome assembly fasta file for each taxon.  There are a number of other parameters of ``ncbi-genome-download`` you can investigate.

    .. code-block:: bash

        for line in `cat ncbi_id.csv`;
            do elements=(${(s:,:)line});
            ncbi-genome-download -s genbank -T ${elements[2]} --verbose --format "fasta,assembly-report" --output ${elements[1]} vertebrate_other; 
        done