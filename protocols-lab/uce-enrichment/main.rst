.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)


UCE Enrichment
==============

:Author: Brant Faircloth, Travis Glenn, John McCormack, Mike Harvey, Laurie Sorenson, Michael Alfaro, Noor White, SI 2012 Seqcap Training Workshop participants, MYcroarray / Arbor Biosciences, Jessie Salter, Carl Oliveros
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.


Modification History
--------------------

.. program-output:: git log --graph --oneline --decorate -- protocols-lab/fix-incorrect-demultiplexing/main.rst


Previous Versions
------------------

+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Version | Description                                                                                                                                                                                                                                                      |
+=========+==================================================================================================================================================================================================================================================================+
| 1.5     | Added on-bead PCR changes for post-enrichment amplification                                                                                                                                                                                                      |
+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1.4     | Update Nextera blockers to give full-length blocking sequence with Inosines as universal blocker.  Recommend final AMPure cleanup at 1.0X versus 1.8X.  This produces larger (on average) contigs following assembly  (March 28, 2013)                           |
+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1.3     | Change probe concentration to 2X of what we originally used.  Update block mix with statement on using custom Cot-1 (e.g. chicken when working with birds).  Title blocker section for TruSeq adapters.  Add an Illumina Nextera blocker section (Sep. 11, 2012) |
+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1.2     | Cleanup to standardize and clarify (Mar. 6, 2012)                                                                                                                                                                                                                |
+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1.1     | Change blocking adapters => 2 pairs of TruSeq primers versus 4.  The new blockers incorporate inosine to bind to a 10 nt index sequence (Nov. 19, 2011)                                                                                                          |
+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1.0     | Original                                                                                                                                                                                                                                                         |
+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Purpose
-------

The purpose of this protocol is to adapt the SureSelect/MYcroarray/Arbor Biosciences enrichment kits for enrichments of UCE DNA libraries prepared using TruSeq/TruSeq-style/Illumina Nextera adapters and common library preparation kits.  We also provide the reagent mixtures used in these kits (from Gnirke et al. 2009 and Blumenstiel et al. 2010), so that you can make more if necessary.

In essence, the protocol provided below is a hybrid of the original UCE enrichment protocol and the protocol provided with version 3 of the MYcroarray / Arbor Biosciences enrichment manual.  Basically, the protocol below is a less stringent version of the MYcroarray protocol.  We have also fleshed out the protocol description with some helpful hits that we've discovered after performing a number of enrichments.

Equipment
---------

* Centrifuge
* Magnet stand or magnet plate
* Water bath or incubator

Materials
---------

* DNA libraries at ~147 ng/uL
* SureSelect or MySelect or IDT bait library (stored at -80 C)
* SureSelect or MySelect hybridization reagents (Box 1 [-20 C] and Box 2 [Room Temperature]) or equivalent hybridization solutions (see below)
* 500 uM adapter blocking mix (IDT DNA) (AKA Block #3; see below)
* Strip tubes and caps or plates and rubber mats
* AMPure XP or Serapure substitute (home-brew AMPure)
* Life Technologies Dynabeads MyOne Streptavidin C1 (Life Technologies 65001)

If you need to prepare additional hybridization reagents:

* 20 X SSPE  (Life Technologies AM9767) (Hyb #1) 
* 0.5 M EDTA (Life Technologies AM9261) (Hyb #2)
* 50 X Denhardt’s Solution (Life Technologies 750018) (Hyb #3)
* 10 % SDS (Life Technologies AM9822) (Hyb #4)
* Human Cot-1 DNA (Life Technologies 15279-101) (Block #1) or Hyblock Cot-1 DNA (`Applied Genetics Laboratories`_)
* Salmon sperm (Life Technologies 15632-011) (Block #2)
* Superase-IN (Life Technologies AM2694) (RNAse Blocker)
* 20 X SSC (Life Technologies AM9770) (for preparing wash buffer)
* 5 M NaCl (Amresco E529-500) (for preparing binding buffer)
	

Preparation
-----------

Pooling Libraries Before Enrichment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on your intended targets/purpose, you may wish to pool several libraries together before enriching the pool of libraries.  This usually makes the process much easier and faster, and we have used this approach extensively when enriching UCE_ libraries for various organismal groups.  We generally pool 8 libraries together at equimolar ratios prior to beginning the enrichment process.

.. figure:: _static/pooling.pdf

   caption



Adapter Blocking Mix
^^^^^^^^^^^^^^^^^^^^

The blocking mix is one of the most critical components during the enrichment.  Without it, you run the risk of adapter-ligated DNA hybridizing together end-to-end ("daisy-chaining").  You pull out what you want, but lots of other stuff you don’t want comes along in a big daisy chain.  Thus, the purpose of Adapter Blocking Mix is to hybridize to the ends of adapter ligated DNA before you add your probe mix.  As such, the blocking mix should match the adapters you’ve added to your libraries.  Over time, we have used several different types of blocking mix, each of which are provided below.

.. attention:: If you are not sure which blocking oligos to use, you need to determine how your libraries were prepared.  Currently, most libraries are dual indexed for illumina, in which case the **8 nucleotide iTru / TruSeq (double indexed)** blockers will work well.


10 nucleotide TruSeq (single indexed) library blocker 
"""""""""""""""""""""""""""""""""""""""""""""""""""""

If you are working with standard Illumina barcodes, you kit may contain the correct blockers for the sequencing adapters.  If you are using longer indexes (e.g. 10 nt), you will need custom blockers.  The blockers below assume 10 nt indexes.  Adjust the number of Inosines (I) to reflect your index length.

1.	You need the following oligos (250 nM synthesis):

.. code-block:: text

	5’- AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT - 3’
	5’- CAAGCAGAAGACGGCATACGAGATIIIIIIIIIIGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT - 3’

2.	Hydrate the above with ddH2O or TLE to 1000 uM (1 times the number of nMol)
3.	Combine 50 uL of each blocker in a 1.5 mL tube
4.	This is now equivalent to Block #3 and contains 500 µM each blocker

8 nucleotide Nextera library blocker 
""""""""""""""""""""""""""""""""""""

If you are enriching libraries prepared using either of Illumina’s Nextera Kits (Illumina Nextera or Illumina Nextera XT), then you need to block indexes on both ends of the library fragments.  Nextera indexes are 8 bp long, each.  The blockers below assume 8 nt indexes (the standard Nextera index length).

1.	You need the following oligos (250 nM synthesis):

.. code-block:: text

    5’ – AATGATACGGCGACCACCGAGATCTACACIIIIIIIITCGTCGGCAGCGTCAGATGTGTATAAGAGACAG - 3’ 
    5’ – CAAGCAGAAGACGGCATACGAGATIIIIIIIIGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG - 3’

2.	Hydrate the above with ddH2O or TLE to 1000 uM (1 times the number of nMol)
3.	Combine 50 uL of each blocker in a 1.5 mL tube
4.	This is now equivalent to Block #3 and contains 500 µM each blocker

8 nucleotide iTru / TruSeq (double indexed) library blocker
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If you are enriching libraries prepared using the iTru indexing system or Illumina's double-indexed TruSeq protocol, you can use the following blocking oligos, which are outward facing and should ensure you do not accidentally produce PCR product from them.

1.	You need the following oligos (250 nM synthesis):

.. code-block:: text

    5’ – AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIIATCTCGTATGCCGTCTTCTGCTTG - 3’ 
    5’ – GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTIIIIIIIIGTGTAGATCTCGGTGGTCGCCGTATCAT - 3’

2.	Hydrate the above with ddH2O or TLE to 1000 uM (1 times the number of nMol)
3.	Combine 50 uL of each blocker in a 1.5 mL tube
4.	This is now equivalent to Block #3 and contains 500 µM each blocker
  	
Bait Mixes
^^^^^^^^^^

When ordering a "custom" enrichment kit, the number of baits synthesized may be greater than the number of baits designed for the targets.  In this case, you can dilute the resulting kit to enrich more samples than the standard amount in a given kit.  Some examples are given below.  However, your mileage may vary, and diluting baits in a custom kit that you are testing for the first time may have negative effects.  It is always best to titrate the baits to attempt to reach some sort of optimal solution.

.. warning:: If you are using a "catalog" kit from MYcroarray / Arbor Biosciences, then you generally **do not want to dilute the bait set** that you have ordered.  These "catalog" kits are already diluted, so they can be sold at a lower price.

Agilent SureSelect
""""""""""""""""""
Depending on how you ordered your probes, they can vary in concentration.  For birds, we have ordered the 2,200 probe set printed 25 times per "array" to fill out the array (55,000 spots / 2,200 probes).  This means that every 1X aliquot of SureSelect probe mix contains 25X the probes that we actually want.  When we target UCEs, we’ve discovered it’s best to use 2X probe mix with each sample or pool of samples.  This is an increase from what we were originally using and is colloquially known as “2X” probe mix.

To enrich 1 sample
++++++++++++++++++

#.	Remove 0.5 uL of MySelect probe mix
#.	Add this to 4.5 uL of RNase free ddH2O

To enrich 12 samples
++++++++++++++++++++

#.	Remove 6.0 uL of MySelect probe mix
#.	Add this to 54.0 uL of RNase free ddH2O

Mycroarray / Arbor Biosciences MYSelect
"""""""""""""""""""""""""""""""""""""""

Similarly, if you’ve printed ~5,000 probes to fill out a 55,000 probe array, then each probe is represented approximately 10 times, so you have a 10X solution of probe mix.  For a single sample, you’ll then want to dilute each aliquot of probes by a factor of 5.

To enrich 1 sample
++++++++++++++++++

1.	Remove 1.0 uL of MySelect probe mix
2.	Add this to 4.0 uL of RNase free ddH2O

To enrich 12 samples
++++++++++++++++++++

1.	Remove 12.0 uL of MySelect probe mix
2.	Add this to 48.0 uL of RNase free ddH2O


Buffers
^^^^^^^

Binding buffer
""""""""""""""

* Assemble the following components in a 50 mL sterile conical tube (note units):

+-------------------------+---------+
| Reagent                 | Amount  |
+=========================+=========+
| 5.0 M NaCl              | 10 mL   |
+-------------------------+---------+
| 1.0 M Tris-HCl (pH 7.5) | 500 µL  |
+-------------------------+---------+
| 0.5 M EDTA              | 100 µL  |
+-------------------------+---------+
| ddH2O                   | 39.4 mL |
+-------------------------+---------+
| Total Volume            | 50.0 mL |
+-------------------------+---------+

* Rotate several times to mix.


Wash Buffer #1
""""""""""""""

* Assemble the following components in a 50 mL sterile conical tube (note units):

+--------------+--------+
| Reagent      | Amount |
+==============+========+
| 20X SSC      | 2.5 mL |
+--------------+--------+
| 10% SDS      | 0.5 mL |
+--------------+--------+
| ddH2O        | 47 mL  |
+--------------+--------+
| Total Volume | 50 mL  |
+--------------+--------+

* Rotate several times to mix.

Wash Buffer #2
""""""""""""""

* Assemble the following components in a 50 mL sterile conical tube (note units):

+--------------+----------+
| Reagent      | Amount   |
+==============+==========+
| 20X SSC      | 250 µL   |
+--------------+----------+
| 10% SDS      | 500 µL   |
+--------------+----------+
| ddH2O        | 49.25 mL |
+--------------+----------+
| Total Volume | 50 mL    |
+--------------+----------+

* Rotate several times to mix.


Steps
-----

Day 1
^^^^^

Toe Pad Wash (for large chunks; can skip small/flaky samples)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#.  Preheat Buffer ATL on heat block to dissolve precipitate
#.  Add 500 ul 100% ethanol to each tube
#.  Incubate samples on a thermomixer at room temperature at 1000 RPM for 5 minutes
#.  Remove ethanol and discard
#.  Add 500 ul 100% ethanol to each tube
#.  Incubate samples on a thermomixer at room temperature at 1000 RPM for 5 minutes
#.  Add 500ul of 1X STE Buffer to each tube
#.  Incubate on a thermomixer at room temperature at 1000rpm for 3-4 hours
#.  Remove STE Buffer and discard
#.  Add 500ul of 1X STE Buffer to each tube
#.  Incubate on a thermomixer at room temperature at 1000rpm for 3-4 hours
#.  Remove STE Buffer and discard
#.  Proceed to **Toe Pad Digestion**


Toe Pad Digestion
"""""""""""""""""

#.  Remove STE Buffer and discard; add the following to each tube containing a toe pad:

        a.  180 ul Buffer ATL
        #.  20 ul Proteinase K

#.  For large chunks, use flame-sterilized forceps to break up toepad as much as possible
#.  Vortex and place in thermomixer at 56˚C at 1000 rpm for ~2 hours
#.  Remove from thermomixer
#.  Using a separate mini pestle for each sample, mash tissue in tubes (REPEAT this step every few hours if necessary – we try and avoid using mini pestles if possible, so as not to lose material)
#.  Vortex and return to thermomixer and incubate at 56˚C at 1000 rpm overnight


Day 2
^^^^^

Toe Pad Digestion
"""""""""""""""""

1.  Remove from thermomixer and add 25 ul 1M DTT to each sample
#.  Vortex, return to thermomixer, and incubate at 56˚C at 1000 rpm for at least 1 hour
#.  Remove from thermomixer and add 15ul proteinase K
#.  Vortex well and place in thermomixer at 56˚C at 1000 rpm for 30 minutes
#.  If tissue is completely digested, move to Step X.  If tissue not digested, continue incubating and smush with mini pestle every few hours until sample is completely digested

Phenol-Chloroform addition
^^^^^^^^^^^^^^^^^^^^^^^^^^

6.  Spin down Phase Lock Gel Light tubes at 12,000 RPM for 30 seconds
#.  Vortex sample after removing from incubation.  Spin down quickly.
#.  Transfer sample to pre-spun Phase Lock Gel tube
#.  In fume hood, add 225ul Phenol:Chloroform:Isoamyl Alcohol (24:25:1)
#.  In fume hood, mix thoroughly by manually rotating tube for 10 minutes
#.  In fume hood, open lid to each tube to vent gas that has built up in each tube
#.  In fume hood, centrifuge at 14,000 RPM for 15 minutes
#.  While tubes are spinning, label a batch of 1.5 mL tubes for final storage (include initial tube number)
#.  In fume hood, pipet the supernatant to final storage tubes while being careful not to disturb the interface between the two layers.  It may be easiest to pour the supernatant into the final tubes rather than pipetting (if you do puncture the interface, return all liquid to the phase lock tube and repeat step 7).  Dispose of Phase Lock Gel Light tubes appropriately.
#.  Add 20 ul 3M NaOAc to each final storage tube
#.  Add 500 ul cold 100% ethanol
#.  Mix well and store in -20˚C for at least 30 minutes.  Tubes can remain at -20˚C overnight, which may increase DNA yields while also, possibly, increasing impurities in the final extract.

Precipitation
"""""""""""""

19.  Remove tubes from -20˚C, and centrifuge at 14,000 RPM for 10 minutes
#.  Pour off supernatant, being careful not to dislodge the DNA pellet
#.  Add 500ul cold 70% ethanol to each sample without disturbing the DNA pellet
#.  Centrifuge at 14,000rpm for 10 minutes
#.  Pour off supernatant and discard, being careful not to dislodge the DNA pellet
#.  Leave the tubes open and dry the DNA pellet in the fume hood (2-3 hours)
#.  Add 50ul 10mM Tris-HCl pH 7.5 to each tube and close the tube.
#.  Store tubes at 4˚C for 24 hrs or bench top overnight
#.  Proceed to DNA quantification
