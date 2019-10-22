.. methuselah documentation master file, created by
   sphinx-quickstart on Thu Oct  3 17:01:11 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Methuselah's documentation!
======================================

Overview:
=========
M
E
T
H 
U ser
S earch for
E pigenetic
L ocus
A ssessment of
H eterogeneity

Optimizing and predicting performance of DNA methylation biomarkers using sequence methylation density information.

Download:
=========
A copy of the package can be obtained by downloading the git repository::

    wget https://github.com/bmill3r/methuselah.git


Installation:
=============
Required dependencies can be found in environment.yml.
Importantly, methuselah is meant for python 3.

samtools 1.6 was used in its development. It is required for running `methuselah READtoMD` on .bam sequence alignment files.


We first recommend installing in a fresh python virtual environment, either mediated by conda or vitrualenv, using the environment.yml.

For conda::

    conda env create -f environment.yml

    conda activate methuselah

Then install::

    cd methuselah

    python setup.py build

    python setup.py install


Check that methuselah is installed with::

    methuselah -V

or::

    methuselah -h



Quick Usage:
============
Generate methylation density (.DT.timestamp.csv) file
-----------------------------------------------------
From either sequencing alignment reads or DREAMing methylation melt data.

For sequences:
^^^^^^^^^^^^^^
i. Each file must be its own .sam or .bam file. Point to the directory that contains the files to analyze.
ii. Alignment files must be from Bismark (https://www.bioinformatics.babraham.ac.uk/projects/bismark/).
iii. Command to process the sequencing reads::

	methuselah READtoMD -i path/to/files --interval BED/with/genome/coordinates

The interval file is a bed file that has genomic coordinates for which methylation density information will be extracted for the samples. If multiple intervals passed, the analysis will be done using pooled information from all of them. Otherwise, can pass "chr:start-stop" for a single genomic region.

For DREAMing data:
^^^^^^^^^^^^^^^^^^
i. Use a raw melt temp .csv file as input. A custom template is provided.
ii. Command to process the DREAMing melt data::

	methuselah DREAMtoMD -i rawmelt.csv --tempsToMDs.csv --numberCpGs

tempsToMDs.csv should be a two column file indicating the methylation density a melting temperature corresponds to for the given locus for which DREAMing was performed. numberCpGs indicates the number of internal CpGs within the given locus for which DREAMing was performed. Both of these are required.

Methylation Density Binary Classifier performance
-------------------------------------------------
Compute the optimal methylation density cutoff for the case and control samples, given the sequence methylation density profiles for the genomic region of interest (or pooled genomic regions::

	methuselah MDBC -i densityTable.csv -a cases -b controls

By default, returns summary tables containing the TPR, 1 - FPR, AUC, and optimal read cutoff for each methylation density cutoff (MDC). Will also return ROCs and boxplots for the optimal MD.

This is done using either the sample read counts or the sample read fractions. As in, the number of reads with a given methylation density or higher, or the sample fraction of reads with a given methylation density or higher.

More extensive plots and sample information can be obtained with additional command flags.

Samples can be normalized based on relative sample fractions, and specific methylation density cutoffs can also be used.

For more information:
---------------------
use::

    methuselah -h

    methuselah READtoMD -h

    methuselah DREAMtoMD -h

    methuselah MDBC -h


Troubleshooting:
^^^^^^^^^^^^^^^^
Overall:
1. Using full, absolute paths is probably a good idea...

Even more detailed information:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Follow the `Documentation` link under 'Contents" below:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Documentation
   License
   Need Help


Indices and tables
==================
* :ref:`search`
