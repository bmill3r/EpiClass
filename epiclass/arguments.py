'''
EpiClass

#########################\\>==+,,++=|\\################################
#######################\,......___,__.-.+\\############################
######################\+,__..___..._+=>>>==\\##########################
#####################\=,____,,++,,_...,=|\\<+>\########################
####################|=,,+++++,,,,+,,,,,,+=<\\\<\#######################
###################|=,,,,,,,,_____,,_,,,,,,>|\\=+<#####################
##################>|+,__...____......._,,,,+|\\|<\#####################
#################\=|+,_..__________....__,,+<\\\\\#####################
################\\>|+,__________.....______+<\\\\\#####################
###############<,_,\|+=====+,___,,,,,,____,<\\\\\\#####################
###############<+_+\|><\\\\\\\\\<\\\||\\||\\\\=,__\####################
###############\=,>\+,,++=<|\=<\+=<<===>\<+||+,..,>####################
################\++<=,,++++\=,_==,++++++=,,>++=__+\####################
#################\,__,,,++=,____,<>++++__,=+_++,+\#####################
##################\_.._..++,_.._,,=,...._,+_.__+#######################
###################\\+___,_________,...__,,_._|########################
#####################\,__,,,,,,,,,___.._,,=\|\#########################
######################|_____,,,_____..._,,<############################
#######################=_____________.___,|############################
########################=,,________,,____+\############################
########################\>,++,,,+,,,,___+\\<\##########################
########################\\>,,+,==___,,=\\\\\+>\\#######################
#####################\\\\\\|++,_,,_+>\\\\\\\>==,_,=|\\#################
################\\\<+<|\\\\\\\>+,++\#\\\\\\\|<<>++,,,,,,+=>|\\\########
##########\\<=+=|\<+<>\###\\\\\\\\\\\\\\\\\\\\|<>>>>>==+++++,++=<\#####

Optimizing and predicting performance of DNA methylation biomarkers using sequence methylation density information.

2019  Brendan F. Miller
bmille79 <at> jh <dot> edu

-----------
  PUBLIC DOMAIN NOTICE
 
  This software is "United States Government Work" under the terms of the United
  States Copyright Act. It was written as part of the authors' official duties
  for the United States Government and thus cannot be copyrighted. This software
  is freely available to the public for use without a copyright
  notice. Restrictions cannot be placed on its present or future use.
 
  Although all reasonable efforts have been taken to ensure the accuracy and
  reliability of the software and associated data, the National Human Genome
  Research Institute (NHGRI), National Institutes of Health (NIH) and the
  U.S. Government do not and cannot warrant the performance or results that may
  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
  disclaim all warranties as to performance, merchantability or fitness for any
  particular purpose.
 
  Please cite the authors in any work or product based on this material.
-------

'''

import sys
import argparse
from textwrap import dedent

from .__init__ import __version__

description_table = """\

    EpiClass

    #########################\\>==+,,++=|\\################################
    #######################\,......___,__.-.+\\############################
    ######################\+,__..___..._+=>>>==\\##########################
    #####################\=,____,,++,,_...,=|\\<+>\########################
    ####################|=,,+++++,,,,+,,,,,,+=<\\\<\#######################
    ###################|=,,,,,,,,_____,,_,,,,,,>|\\=+<#####################
    ##################>|+,__...____......._,,,,+|\\|<\#####################
    #################\=|+,_..__________....__,,+<\\\\\#####################
    ################\\>|+,__________.....______+<\\\\\#####################
    ###############<,_,\|+=====+,___,,,,,,____,<\\\\\\#####################
    ###############<+_+\|><\\\\\\\\\<\\\||\\||\\\\=,__\####################
    ###############\=,>\+,,++=<|\=<\+=<<===>\<+||+,..,>####################
    ################\++<=,,++++\=,_==,++++++=,,>++=__+\####################
    #################\,__,,,++=,____,<>++++__,=+_++,+\#####################
    ##################\_.._..++,_.._,,=,...._,+_.__+#######################
    ###################\\+___,_________,...__,,_._|########################
    #####################\,__,,,,,,,,,___.._,,=\|\#########################
    ######################|_____,,,_____..._,,<############################
    #######################=_____________.___,|############################
    ########################=,,________,,____+\############################
    ########################\>,++,,,+,,,,___+\\<\##########################
    ########################\\>,,+,==___,,=\\\\\+>\\#######################
    #####################\\\\\\|++,_,,_+>\\\\\\\>==,_,=|\\#################
    ################\\\<+<|\\\\\\\>+,++\#\\\\\\\|<<>++,,,,,,+=>|\\\########
    ##########\\<=+=|\<+<>\###\\\\\\\\\\\\\\\\\\\\|<>>>>>==+++++,++=<\#####

    Optimizing and predicting performance of DNA methylation biomarkers using methylation density information.

    Copyright (C) 2019  Brendan F. Miller
    bmille79 <at> jh <dot> edu

    This program is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with
    this program.  If not, see <https://www.gnu.org/licenses/>.

    ---------------------------------------------------------
    ---------------------------------------------------------
    ---------------------------------------------------------
    
    The epiclass sub-modules can be accessed by executing:
        'epiclass module_name arg1 arg2 ...'
    Sub-module help  and argsuments can be displayed by executing:
    'epiclass module_name --help'
    
    Sub-module descriptions:
    
        +---------------+---------------------------------------------------------------------------------------+
        |   DREAMtoMD   |   Convert table of raw DREAMing melting peak counts to methylation density table      |
        |---------------|---------------------------------------------------------------------------------------|
        |   READtoMD    |   Convert .bam/.sam bismark alignment files to methylation density table              |
        |---------------|---------------------------------------------------------------------------------------|
        |               |   Identify the optimal methylation density (MD) and read value cutoffs between cases  |
        |     MDBC      |   and controls. Returns summary table and plots of optimal MD cutoff. Additional      |
        |               |   visualization and information collection options available.                         |
        +---------------+---------------------------------------------------------------------------------------+
    
    ---------------------------------------------------------
    Notes for DREAMtoMD:
    IMPORTANT
    Converting the data into 'methylationDensityTable.csv' requires a "tempsToMDs" table, and "numCpGs" parameters
    that have already been determined for the DREAMing locus in question.
    
    The input CSV file containing raw DREAMing melting peak temperature data:
    Structure should be as follows:
    sample      XXXX_L    XXXX_H 
    copies_loaded    6000    6000
    1       80.8
    2       80.8
    3       80.8    84.0
    4       80.8
    5       80.8
    6       80
    7       80.8    83.8
    8       80.8    82.6
    9       80.8
    10      80
    11      80
    12      80.8
    plate       A       A
    date        9.12    9.12
    sample --> first row are the samples names. Each sample has two columns designated by '_L' and '_H'. L = lower melting peak temperature and H = temperature of second melting peak farthest to the right (if exists) (aka the highest on the melt trace).
    copies_loaded --> row containing number of genomic equivalents loaded into the DREAMing assay for the given sample on the given plate and date run. Because two columns are designated for a sample, value should be the same for each of those cells.
    1:12 --> wells 1-12 on a given row of a 96-well microtiter plate for which the DREAMing assay was run for the given sample. Note that the number of wells can vary depending on the passay or plate used.
    plate --> optional row assigning samples to the DREAMing assay for which they were run. Because two columns are designated for a sample, value should be the same for each of those cells.
    date --> optional row assigning samples to the date for which they were run in a given DREAMing assay. Because two columns are designated for a sample, value should be the same for each of those cells.
    Additional rows could be added.
    Note that a given sample can have several replicates, in which the column name will then be followed by a .1, .2, etc depending on the number of replicates. ex: XXXX_L.1, XXXX_H.1; XXXX_L.2, XXXX_H.2
    By default, replicates of a sample are combined to make the 'plasmaMeltPeakCounts.csv' or 'methylationDensityTable.csv' files.
    IMPORTANT:
    fractionless melting peak temperatures (80.0C, ex.) should be recorded as '80' in order to match 'meltTempsToMD.csv'
    
    tempsToMDs:
    ex:
    For ZNF154 DREAMing assay:
    '79.2':0.0,
    '79.4':0.0,
    '79.6':0.0,
    '79.8': 0.0,
    '80': 0.07,
    '80.2': 0.14,
    '80.4': 0.14,
    '80.6': 0.21,
    '80.8': 0.28,
    '81': 0.35,
    '81.2': 0.35,
    '81.4': 0.42,
    '81.6': 0.49,
    '81.8': 0.56,
    '82': 0.56,
    '82.2': 0.63,
    '82.4': 0.7,
    '82.6': 0.77,
    '82.8': 0.84,
    '83': 0.84,
    '83.2': 0.91,
    '83.4': 1.0,
    '83.6': 1.0,
    '83.8': 1.0,
    '84': 1.0
    
    Input2Bg = False ignores the background input copies used for loading estimating number of Genomic DNA Equivalents loaded into a DREAMing assay.
    As such, EFs could be fractions of reads relative to the recorded reads in a sample and not adjusted to the estimated Genomic DNA Equivalents loaded.
    
    ---------------------------------------------------------
    Notes for READtoMD:
    Alignment files are sorted alignment files from Bismark. If BAM, requires samtools view to access the reads.
    Otherwiise, will need alignment files to be in SAM.
    Extracted Methylation sequences come from reads that overlap intervals in interest. Sequence can be entire sequence or only sequence in overlap from a read or stitched together from paired-reads.
    Example of a line (read) from a bismark alignment file that is expected:
    readID 163 chrom readStart 40 120M = 74960987 158 actualDNASequence alignmentScore NM:i:30 MD:Z: XM:Z:methylationInfo(...z...h...x..Z...etc) XR:Z:GA XG:Z:GA XB:Z:6
    Methylation information taken from the 'XM:Z:' line position.
    Requires a directory that contains the alignment files of interest. Each SAM/BAM file will be considered a separate sample to include in the methylationDensityTable.csv.
    Requires a tab-delimited BED file where each row contains the chromosome, start, and stop coordinates of interval(s) of interest to extract reads from. (no headers)
    Can also be a single command line input (chr:start-stop)

    methylationDensityTable.csv columns:
    ['chr', 'interval_start', 'interval_stop', 'methSeqStart', 'methSeqStop', 'methSeqLength', 'numU', 'numM', 'MD', sample.1, sample.2,...]
    chr = chromosome location of reads in file
    interval_start = start of a interval of interest that contains the methylation sequence.
    interval_stop = stop of a interval of interest that contains the methylation sequence.
    methSeqStart = start position of the methylation sequence extracted.
    methSeqStop = stop position of the methylation sequence extracted.
    methSeqLength = length of the methylation sequence of interest.
    numU = number of unmethylated CpGs determined by bismark alignment of the extracted methylation sequence.
    numM = number of methylated CpGs determined by bismark alignment of the extracted methylation sequence.
    MD = 'methylation density'; numM / (numU + numM)
    sample = column with counts of methylation sequences in sample that have the specific information supplied in previous columns.
    
    ---------------------------------------------------------
    Notes for MDBC:
    Returns summary table of AUC, optimal read value cutoff, sensitivity, and specificity for each MD cutoff.
    Read value can be normalized sample read counts or normalized sample read fractions for the MD in question.
    Returns ROC curve and boxplots of overall optimal MD cutoff for cases vs controls.
    If indicated, can return additional plots or average methylation analysis. See flags.
    Normalized first to sample input fractions, if given. Could be number of arbitrary units of sample analyed or fraction of sample relative to other samples.
    Then, read counts or fractions normalized to number of CpGs covered by each read.

    Some calculations defined:
    Average methylation:
    number of meCpGs in sample reads / total CpGs in all sample reads
    EF = epiallelic fraction. Defined as:
    number of CpGs covered by sample reads with MD >= cutoff / total CpGs in all sample reads
    (special case: if MD cutoff = 0, then choose reads with MD > 0)
    
    Table with relative sample input fractions.
    If supplied, then will perform sample adjustments based on fractions in table.
    First column is sample names ('samples') of a given sample set and second column is sample fractions ('fractions') within the set.
    ex:
    samples    fractions
    A       1.0
    B       0.75
    C       0.95
    ...
    OR
    samples    fractions
    A       55.0
    B       102.6
    C       20.0
    ...
    where the fractions now represent actual mg amounts of stool (for example) processed per sample.
    The fractions could all be based on the relative fraction of starting material used.
    As in, if 500mg of stool was used for each sample, then the fractions could be the fraction of this
    500mg used per sample. Thus sample read counts for fractions would be normalized based on 500mg of stool.
    Alternatively, fractions could be changed to the actual amount of material used.
    For example, 50mg stool, 125mg stool, 75.5mg of stool; or 1mL plasma, 0.75mL plasma, etc.
    Then, read counts or fractions would be reads/sampleInput, where sampleInput could be mg stool, mLs plasma, etc.
"""

epilog_description = '''\

    example usage:
    
    ./epiclass DREAMtoMD -i ./rawDREAMing.csv --tempsToMDs tempsToMDs.csv --numberCpGs 14 -o ./label_
    
    will return ./label_rawDREAMing.MDT.csv methylation density table. Then, input this into:
    
    ./epiclass MDBC -i -a cases -b controls --fractions ./sampleInputFractions.csv
    
    which will return:
    ./MDBC.label_rawDREAMing.READ-COUNT-SUMMARY.csv, ./MDBC.label_rawDREAMing.EF-SUMMARY.csv,
    ./MDBC.label_rawDREAMing.READ-COUNT-OPTMD-BOX.png, ./MDBC.label_rawDREAMing.EF-OPTMD-BOX.png,
    ./MDBC.label_rawDREAMing.READ-COUNT-OPTMD-ROC.png, ./MDBC.label_rawDREAMing.EF-OPTMD-ROC.png
    
    where *summary.csv files are tables of TPR, 1-FPR, AUC, optimal read value cutoffs for each MD cutoff assesed
    and *BOX.png and *ROC.png are boxplots or ROC curves using the sample read values for the optimal MD cutoff.

'''

def get_arguments(args):

    if args is None:
        args = sys.argv[1:]

    """INITIALIZE PARSER"""
    parser = argparse.ArgumentParser(
        prog='epiclass',
        description=dedent(description_table),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=dedent(epilog_description))
    # Optional arguments
    parser.add_argument(
        '-V', '--version',
        help='Print installed version to stout',
        action='version',
        version='%(prog)s ' + str(__version__))

    """MODULE SUBPARSER PROGRAMS"""
    subparser = parser.add_subparsers(dest='cmd')

    """DREAMtoMD SUBPARSER"""
    dreamer = subparser.add_parser(
        'DREAMtoMD',
        description='convert raw DREAMing data table to methylation denstiy table',
        add_help=False)
    # Required arguments
    dreamer_reqs = dreamer.add_argument_group('required arguments')
    dreamer_reqs.add_argument(
        '-i', '--input',
        help='Path to raw DREAMing table (FULL path)',
        metavar='<path>',
        type=str,
        required=True)
    dreamer_reqs.add_argument(
        '-temps', '--tempsToMDs',
        help='''REQUIRED.
        Path to table with melting peak temperature to assigned methylation density value for the given DREAMing locus.
        Used to calibrate the melting peak temperature to methylation density.
        Generated by modeling the observed melting temperature and the actual number of meCpGs via amplicon sequencing.
        IMPORTANT: the temperatures in the raw DREAMing data file must fall within the range of melting temperatures in this file.''',
        metavar='/meltTempsToMD.csv',
        type=str,
        required=True)
    dreamer_reqs.add_argument(
        '-cpg', '--numberCpGs',
        help='''REQUIRED.
        Number of CpGs covered in DREAMing locus and previously used to calibrate the DREAMing melting peak temperature to methylation density.''',
        metavar='#CpGs',
        type=int,
        required=True)
    # Optional arguments
    dreamer_opts = dreamer.add_argument_group('optional arguments')
    dreamer_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    dreamer_opts.add_argument(
        '-v', '--verbose',
        help='Set verbosity level. -v=WARNING, -vv=INFO, -vvv=DEBUG',
        action='count',
        default=0,
        required=False)
    dreamer_opts.add_argument(
        '--input2bg',
        help='Include input genomic equivalents as background reads with a methylation density = 0.',
        action='store_true',
        default=False,
        required=False)
    dreamer_opts.add_argument(
        '--poisson',
        help='Adjust the counts of melting peak temperatures for each sample in a given DREAMing assay based on a Poissonian distribution.',
        action='store_true',
        default=False,
        required=False)
    dreamer_opts.add_argument(
        '--tempResolution',
        help='''Melting temperature resolution of the melting peaks.
        Typically 0.2C, but could change depending on the thermocycler used for the DREAMing assay.''',
        metavar='0.2',
        type=float,
        default=0.2,
        required=False)
    dreamer_opts.add_argument(
        '--includedRows',
        help='''Comma separated list of rows of interest for the analysis.
        'Sample', and 'copies_loaded' required.
        However, additional labels could be added to include other sample/assay information in the file output table.''',
        metavar='Sample,copies_loaded',
        type=str,
        default='Sample,copies_loaded',
        required=False)
    dreamer_opts.add_argument(
        '-o', '--output',
        help='Path to output directory',
        metavar='<path>',
        type=str,
        required=False)

    """READtoMD SUBPARSER"""
    reader = subparser.add_parser(
        'READtoMD',
        description='convert sam/bam bismark aligned reads to methylation denstiy table',
        add_help=False)
    # Required arguments
    reader_reqs = reader.add_argument_group('required arguments')
    reader_reqs.add_argument(
        '-i', '--input',
        help='Path to directory containing bam/sam alignment files (one per sample of interest) (FULL path)',
        metavar='<path>',
        type=str,
        required=True)
    reader_reqs.add_argument(
        '--intervals',
        help='''Path to tab-delimited 'intervalFile.txt' or BED file (no column names) genomic regions.
        Or 'chr:start-stop'. ''',
        metavar='''<path> or chr:start-stop''',
        type=str,
        required=True)
    # Optional arguments
    reader_opts = reader.add_argument_group('optional arguments')
    reader_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    reader_opts.add_argument(
        '-v', '--verbose',
        help='Set verbosity level. -v=WARNING, -vv=INFO, -vvv=DEBUG',
        action='count',
        default=0,
        required=False)
    reader_opts.add_argument(
        '--fileType',
        help='''Label indicating the type of files to select.
        Typically contains information like the locus.
        For example: 'ZNF154reads.sam' will indicate selection of the *ZNF154reads.sam* files in the input directory.
        i.g., HCC1.dedup.SOcoord.ZNF154reads.sam (sample = HCC1)''',
        metavar='label',
        type=str,
        required=False)
    reader_opts.add_argument(
        '--overlap',
        help='''Extract methylation sequences from read/read pairs that overlap the defined interval at least this much.''',
        metavar='0.0',
        default=0.0,
        type=float,
        required=False)
    reader_opts.add_argument(
        '--slice',
        help='''Flag call to use methylation sequences information only from sections of selected reads that are contained in the overlap with the interval(s).
        Otherwise use methylation information from the entire read(s) that overlap.''',
        action='store_true',
        default=False,
        required=False)
    reader_opts.add_argument(
        '--fileTag',
        help='''Optional label to tag output files with. Otherwise a timestamp is given.''',
        metavar='label',
        type=str,
        required=False)
    reader_opts.add_argument(
        '-o', '--output',
        help='Path to output directory',
        metavar='<path>',
        type=str,
        required=False)

    """MDBC SUBPARSER"""
    mdbcer = subparser.add_parser(
        'MDBC',
        description='Methylation Density Binary Classifier (MDBC) analysis with plotting.',
        add_help=False)
    # Required arguments
    mdbcer_reqs = mdbcer.add_argument_group('required arguments')
    mdbcer_reqs.add_argument(
        '-i', '--input',
        help='Path to methylation density table.',
        metavar='<path>',
        type=str,
        required=True)
    mdbcer_reqs.add_argument(
        '-a', '--cases',
        help='''List of names of samples representing the case set.
        Should match subset of sample columns in the input methylation density table.
        Can also enter single string to use as ID for selecting samples in density table as cases.''',
        metavar='case1',
        type=str,
        nargs='+',
        required=True)
    mdbcer_reqs.add_argument(
        '-b', '--controls',
        help='''List of names of samples representing the control set.
        Should match subset of sample columns in the input methylation density table.
        Can also enter single string to use as ID for selecting samples in density table as controls.''',
        metavar='control1',
        type=str,
        nargs='+',
        required=True)
    # Optional arguments
    mdbcer_opts = mdbcer.add_argument_group('optional arguments')
    mdbcer_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    mdbcer_opts.add_argument(
        '-v', '--verbose',
        help='Set verbosity level. -v=WARNING, -vv=INFO, -vvv=DEBUG',
        action='count',
        default=0,
        required=False)
    mdbcer_opts.add_argument(
        '--fractions',
        help='''Path to table with relative sample fractions.''',
        metavar='<path>',
        type=str,
        required=False)
    mdbcer_opts.add_argument(
        '--MDcutoffs',
        help='''List of methylation density cutoff values as fraction out of 1.
        Use to assess the perfomance of the MDBC at each one.''',
        metavar='0.0 0.05',
        type=float,
        nargs='+',
        required=False)
    mdbcer_opts.add_argument(
        '--maxCutoff',
        help='''Optionally set the maximum read count cutoff or EF cutoff range to test.''',
        metavar='100.0',
        type=float,
        required=False)
    mdbcer_opts.add_argument(
        '--hdf_label',
        help='''Path of HDF5 file which stores samples values (normalized read counts or fractions) for each MD cutoff.''',
        metavar='<path>',
        type=str,
        required=False)
    mdbcer_opts.add_argument(
        '--ignoreCountsummary',
        help='Skip processing read count summary performance and optimal MD cutoff.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--ignoreEFsummary',
        help='Skip processing read fraction summary performance and optimal MD cutoff.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--totalReadCounts',
        help='Return table of normalized total methylated read counts for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--totalReadCountPlot',
        help='Return boxplot of normalized total methylated read counts for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--totalEf',
        help='Return table of normalized total methylated read count fractions for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--totalEfPlot',
        help='Return boxplot of normalized total methylated read count fractions for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--readcountsEachMD',
        help='Return table of normalized read counts for each MD cutoff for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--EfEachMD',
        help='Return table of normalized read fractions for each MD cutoff for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--sampleValsAtMD',
        help='Return tables of normalized read fractions and efs for each MD cutoff value indicated in list for cases and controls.',
        metavar='0.0 0.05',
        type=float,
        nargs='+',
        required=False)
    mdbcer_opts.add_argument(
        '--sampleAveMethTable',
        help='Return table of average methylation values for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--casesReadCountPlot',
        help='Return stacked barplot showing normalized read counts for each MD cutoff in cases.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--controlsReadCountPlot',
        help='Return stacked barplot showing normalized read counts for each MD cutoff in controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--casesEfPlot',
        help='Return stacked barplot showing normalized sample read fractions for each MD cutoff in cases.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--controlsEfPlot',
        help='Return stacked barplot showing normalized sample read fractions for each MD cutoff in controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--readcountDistributionPlot',
        help='Return histogram showing distribution of normalized read count proportions for each MD for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--EfDistributionPlot',
        help='Return histogram showing distribution of normalized read fractions proportions for each MD for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--readcountheatmaps',
        help='Return heatmaps of TRP, FPR, and TPR-FPR for each MD cutoff and normalized read count cutoff.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--Efheatmaps',
        help='Return heatmaps of TRP, FPR, and TPR-FPR for each MD cutoff and normalized read fraction cutoff.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--optimalMDreadcounts',
        help='Return values of normalized sample read counts for reads with methylation densities at or above the optimal MD cutoff.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--optimalMDEf',
        help='Return values of normalized sample read fractions for reads with methylation densities at or above the optimal MD cutoff.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--fileTag',
        help='''Optional label to tag output files with.''',
        metavar='label',
        type=str,
        required=False)
    mdbcer_opts.add_argument(
        '-o', '--output',
        help='Path to output directory',
        metavar='<path>',
        type=str,
        required=False)

    args = parser.parse_args(args)
    args_dict = vars(args)

    return args, args_dict
