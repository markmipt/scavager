Scavager - a proteomics post-search validation tool
---------------------------------------------------

The pepXML or MzIdentML files are required for basic operation of the script. Currently supported search engines:
IdentiPy, X!Tandem, Comet, MSFragger, MSGF+, Morpheus.

FASTA file is required for calculation of NSAF (label-free quantitation index), protein sequence coverage and amino acid statistics.

For MSGF+ and Morpheus search engines it is desirable to provide cleavage rules used in search
(these search engines do not report number of missed cleavages for peptides).

The output of Scavager contains:

- tab-separated table with unfiltered peptide-spectrum matches (ends with _PSMs_full.tsv)
- tab-separated table with identified peptide-spectrum matches at 1% PSM FDR (ends with _PSMs.tsv)
- tab-separated table with identified peptides at 1% peptide FDR (ends with _peptides.tsv)
- tab-separated table with identified proteins without grouping at 1% protein FDR (ends with _proteins.tsv)
- tab-separated table with identified protein groups at 1% protein FDR (ends with _protein_groups.tsv)
- PNG figure with PSM, peptide and protein features distributions

Citing Scavager
---------------
Ivanov et al. Scavager: A Versatile Postsearch Validation Algorithm for Shotgun Proteomics Based on Gradient Boosting. doi: 10.1002/pmic.201800280

Installation
------------
Using pip:

    pip install Scavager


Usage
-----
Algorithm can be run with following command (works with Python2.7/Python3+):

    scavager path_to_pepXML/MZID

OR

    scavager -h

Protein grouping using DirectMS1 results
----------------------------------------
Protein groups can be generated using parsimony principle combined with information from MS1 spectra:

    scavager path_to_pepXML/MZID -ms1 path_to_DirectMS1_proteins_full_noexclusion.tsv

Details on combination of parsimony principle and MS1 information are available at: https://github.com/markmipt/protein_inference_using_DirectMS1

Protein grouping for indistinguishable proteins
------------------------------------------------
By default, when multiple proteins have the same sets of peptides, the Scavager choose protein group leader using alphabetical order. However, it is possible
to choose group leader randomly by using "-sr" option. The same option can be used with MS1 spectra information if multiple proteins have both same sets of MS/MS identifications and DirectMS1 scores.


Links
-----

- GitHub repo & issue tracker: https://github.com/markmipt/scavager
- Mailing list: pyteomics@googlegroups.com, markmipt@gmail.com
