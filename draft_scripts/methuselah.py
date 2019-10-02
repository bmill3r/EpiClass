#!/usr/bin/env	python

import sys
import os
import os.path
from os import walk
import subprocess
from distutils.spawn import find_executable

import re
from collections import Counter
import itertools
import datetime

import numpy as np
import pandas as pd
import tables

import scipy.stats as stats
from sklearn.metrics import roc_curve, auc
from matplotlib import pyplot as plt

import argparse
from textwrap import dedent

import matplotlib as mpl
#aesthetics for plots:
font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 12}
mpl.rc('font', **font)

now = datetime.datetime.now()
timestamp = now.strftime("%Y-%m-%d_%H-%M")
date = str(datetime.date.today())
cwd = os.getcwd()

# print(std.out 'print' statements to new file and to console:)
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open('methuselah.' + date + '.log', "a+")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        pass

sys.stdout = Logger()

# get file name from path no matter the operating system or file path
def path_leaf(path):
    head, tail = os.path.split(path)
    return tail or os.path.basename(head)

description_table = """\

    The methuselah sub-modules can be accessed by executing:
        'methuselah module_name arg1 arg2 ...'
    Sub-module help  and argsuments can be displayed by executing:
    'methuselah module_name --help'
    
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
    
    ./methuselah DREAMtoMD -i ./rawDREAMing.csv --tempsToMDs tempsToMDs.csv --numberCpGs 14 -o ./label_
    
    will return ./label_rawDREAMing.MDT.csv methylation density table. Then, input this into:
    
    ./methuselah MDBC -i -a cases -b controls --fractions ./sampleInputFractions.csv
    
    which will return:
    ./label_rawDREAMing.MDT.countsummary.csv, ./label_rawDREAMing.MDT.efsummary.csv,
    ./label_rawDREAMing.MDT.readcount.optMD.BOX.png, ./label_rawDREAMing.MDT.ef.optMD.BOX.png,
    ./label_rawDREAMing.MDT.readcount.optMD.ROC.png, ./label_rawDREAMing.MDT.ef.optMD.ROC.png
    
    where *summary.csv files are tables of TPR, 1-FPR, AUC, optimal read value cutoffs for each MD cutoff assesed
    and *BOX.png and *ROC.png are boxplots or ROC curves using the sample read values for the optimal MD cutoff.

'''

__version__ = '2.0.0-beta'


def get_arguments(args):

    if args is None:
        args = sys.argv[1:]

    """INITIALIZE PARSER"""
    parser = argparse.ArgumentParser(
        prog='methuselah',
        description=dedent(description_table),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=dedent(epilog_description))
    # Optional arguments
    parser.add_argument(
        '-v', '--version',
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
        help='Path to raw DREAMing table',
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
        help='Path to directory containing bam/sam alignment files (one per sample of interest)',
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
        '--totalreadcounts',
        help='Return table of normalized total methylated read counts for cases and controls.',
        action='store_true',
        default=False,
        required=False)
    mdbcer_opts.add_argument(
        '--totalreadcountsPlot',
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
        '-o', '--output',
        help='Path to output directory',
        metavar='<path>',
        type=str,
        required=False)

    args = parser.parse_args(args)
    args_dict = vars(args)

    return args, args_dict

class dreamingToDensityTable():

    def __init__(self, rawDreamingData, tempsToMDs, numberCpGs,
                 meltResolution=0.2, includeBackground=False,
                 poissonAdjustment=False, includedRows='Sample,copies_loaded'):

        self.rawDreamingData = pd.read_csv(rawDreamingData, index_col=[0])

        calibrationDF = pd.read_csv(tempsToMDs)
        tempValues = [str(int(i)) if str(i)[-1] == '0' else str(i)
                      for i in calibrationDF['temp']]
        MDvalues = [float(i) for i in calibrationDF['MD']]
        self.tempsToMDs = dict(zip(tempValues, MDvalues))

        # if counting input copies as part of background, then consider as fragments with MD = 0%:
        if includeBackground == True:
            self.tempsToMDs['copies_loaded'] = 0.0

        self.numberCpGs = numberCpGs
        self.numM = np.arange(0, self.numberCpGs+1)
        self.numU = np.arange(0, self.numberCpGs+1)[::-1]

        self.meltResolution = meltResolution
        self.includedRows = list(includedRows.split(','))

        self.poissonAdjustment = poissonAdjustment

        wells = []
        for i in self.rawDreamingData.index.values:
            try:
                wells.append(int(i))
            except:
                ValueError
        self.numWells = len(wells)

        # in table, "copies_loaded" = row 0, and wells start at row 1. Additional info rows are after well rows
        self.wellIndices = [str(i) for i in np.arange(1, self.numWells+1)]

        # list of unique sample names
        # split on _H or _L headings that designate type of melt peak in raw DREAMing table
        self.samplesOfInterest = sorted(
            list(set([re.split(r'_[HL]', i)[0] for i in self.rawDreamingData.columns])))

        # get lowest and highest melting temperatures recored in table to set range
        self.minMelt = np.nanmin(
            self.rawDreamingData.loc[self.wellIndices].astype(float).values.ravel())
        self.maxMelt = np.nanmax(
            self.rawDreamingData.loc[self.wellIndices].astype(float).values.ravel())

        # range of recorded melting peak temperatures in raw DREAMing data to use as temporary column headings
        # to store melting peak temperature counts for samples.
        self.meltCols = [str(round(t, 1)) if str(round(t, 1))[-1] != '0' else str(round(t, 1))[:2]
                         for t in np.arange(self.minMelt, self.maxMelt+self.meltResolution, self.meltResolution)]

        # standard columns of interest
        self.columns = self.includedRows + self.meltCols

        # index for appending melting peak counts and other info for each sample to proper location in sampleRecord
        self.recordIndex = dict(
            zip(self.columns, np.arange(len(self.columns))))

    @property
    def getMeltPeakCounts(self):
        '''
        Make list for each sample where each position corresponds to associated entry in the recordIndex.
        Update this sample list with counts of melting peak temperatures or other sample information.
        Append these sample lists as 'rows' to master dataframe:
        '''
        master_array = []
        for sample in self.samplesOfInterest:

            # make list to update for each sample
            sampleRecord = [0] * len(self.columns)

            # get dataframe slice where columns represent melting peak temperature records for instances of sample
            # could be multiple records of a given sample
            # and each sample will have at least two columns: one for L and one for H melt peak recordings
            temp = self.rawDreamingData[self.rawDreamingData.columns[self.rawDreamingData.columns.str.contains(
                sample)]]

            # for each instance of the sample, get rows with melt peak recordings
            # and get number of occurances of each melt peak
            for col in temp.columns:

                meltPeaks = temp.loc[self.wellIndices, col]
                # replace empty cells with 0
                meltPeaks.fillna(0, inplace=True)
                # convert melt peak values to strings
                meltPeaksToStr = meltPeaks.astype(str)

                types, counts = np.unique(meltPeaksToStr, return_counts=True)

                # for each melting peak temp record and count, check to make sure string (ex: '79.8' or '80')
                # and not nan
                for j in np.arange(len(counts)):

                    # ignore the cells that were empty
                    if types[j] != '0':

                        # append counts to correct position in sampleRecord based on recordIndex

                        # if all wells are that melt temp
                        if int(counts[j]) == self.numWells:
                            sampleRecord[self.recordIndex[types[j]]
                                         ] += self.numWells

                        else:
                            # if poisson adjustment:
                            if self.poissonAdjustment == True:
                                sampleRecord[self.recordIndex[types[j]]] += round(
                                    (-np.log(1.0 - (counts[j]/float(self.numWells))) * float(self.numWells)), 0)
                            # or just append number of counts of the given melt temp:
                            else:
                                sampleRecord[self.recordIndex[types[j]]
                                             ] += counts[j]

            # append other info to sampleRecord:
            for i in self.includedRows:

                if i == 'copies_loaded':
                    # update list with total copies loaded for all DREAMing assays of the sample
                    # copies_loaded records should actually be string values in table
                    sampleCopies = [
                        float(val) for val in temp.loc[i, :].values if type(val) == str]

                    # sum together total copies for all DREAMing assay instances of sample
                    # divide by 2 b/c L and H columns count the copies twice
                    sampleRecord[self.recordIndex[i]] = sum(sampleCopies)/2.0

                if i == 'Sample':
                    # update list with sample name
                    sampleRecord[self.recordIndex[i]] = sample

            # append sample row to the master array:
            master_array.append(sampleRecord)

        master_df = pd.DataFrame(master_array, columns=self.columns)
        return master_df

    def densityTable(self):

        # set sample names as index
        df = self.getMeltPeakCounts.set_index('Sample')
        # set columns as methylation density values
        df = df[self.tempsToMDs.keys()]
        # convert melt temperature column names to MD values
        df.rename(index=str, columns=self.tempsToMDs, inplace=True)
        # then transpose sample names to columns
        methDensTable = df.T
        # combine the counts for replicate MD recordings
        methDensTableCombined = methDensTable.groupby(
            methDensTable.index).sum()
        # Add in numU, numM, and MD columns
        methDensTableCombined['numU'] = self.numU
        methDensTableCombined['numM'] = self.numM

        MDs = methDensTableCombined.index.values.tolist()
        methDensTableCombined['MD'] = MDs

        methDensTableCombined.reset_index(inplace=True, drop=True)
        # reorder columns:
        methDensTableCombined = methDensTableCombined[[
            'numU', 'numM', 'MD'] + self.samplesOfInterest]

        return methDensTableCombined

class readsToDensityTable():

    
    def __init__(self, fileDirectory, intervals, fileType=None, readSlice=False, overlap=0.0):

        # path to folder containing sam/bam read files
        if str(fileDirectory).endswith('/'):
            self.fileDirectory = fileDirectory
        else:
            self.fileDirectory = str(fileDirectory) + '/'

        # selected sequencing read files
        self.files = []
        for (dirpath, dirnames, filenames) in walk(self.fileDirectory):
            filenames = [dirpath + '/' + f for f in filenames]
            self.files.extend(filenames)

        # filter for .bam/.sam extention files
        self.files = [f for f in self.files if bool(
            re.search(r'[sb]am$', path_leaf(f), re.IGNORECASE))]

        # OPTINAL: only process files that contain fileType str label
        if type(fileType) == str:
            self.files = [
                f for f in self.files if fileType in path_leaf(f)]

        # check if interval file used, otherwise use single interval input of chr:start-stop
        if os.path.exists(intervals) == True:
            self.intervalFileName = intervals
            self.intervals = pd.read_table(self.intervalFileName, names=[
                'chr', 'start', 'stop'], index_col=False)
        else:
            chrom = intervals.split(':')[0]
            start = int(intervals.split(':')[1].split('-')[0])
            stop = int(intervals.split(':')[1].split('-')[1])
            self.intervalFileName = str(cwd) + '/' + 'interval.txt'
            self.intervals = pd.DataFrame(
                {'chr': [chrom], 'start': [start], 'stop': [stop]})

            self.intervals.to_csv(self.intervalFileName, columns=[
                'chr', 'start', 'stop'], header=False, index=False, sep='\t')

        self.readSlice = readSlice
        self.overlap = overlap

    def samtoolsReadLines(self, inputFile):

        lines = []
        cmd = 'samtools view -L ' + \
            str(self.intervalFileName) + ' ' + str(inputFile)
        reads = subprocess.Popen(
            [cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
        stdout, stderr = reads.communicate()

        for line in stdout.splitlines():
            lines.append([line.split('\t')[i]
                          for i in [0, 2, 3, 13]])
        return lines

    def manualReadLines(self, inputFile):

        lines = []
        with open(inputFile, 'r') as f:
            mylist = f.read().splitlines()
            for line in mylist:
                # ignore header lines
                if line.startswith('@'):
                    continue
                else:
                    # read_name, chr, read_start, meth
                    lines.append([line.split('\t')[i] for i in [0, 2, 3, 13]])
        return lines

    def overlappingReads(self, read_df, interval):

        # collect methylation information for sample reads that overlap interval
        interval_reads = []

        int_chr = interval.chr
        int_start = interval.start
        int_stop = interval.stop
        intervalname = str(interval.chr) + ':' + \
            str(int(interval.start)) + '-' + str(int(interval.stop))

        # keep reads that are on same chromosome and overlap interval
        readsIn = read_df[(read_df['chr'] == int_chr)]
        readsIn = readsIn[~(int_start >= readsIn['read_end'])]
        readsIn = readsIn[~(readsIn['read_start'] >= int_stop)]
        readsIn['int_start'] = int_start
        readsIn['int_stop'] = int_stop

        # get length and proportion of read overlap to corresponding interal
        Ecols = readsIn[['read_end', 'int_stop']].idxmin(axis=1)  # overlap end
        readsIn['E'] = [readsIn.iloc[i][Ecols.iloc[i]]
                        for i in np.arange(0, len(Ecols))]
        Scols = readsIn[['read_start', 'int_start']].idxmax(
            axis=1)  # overlap start
        readsIn['S'] = [readsIn.iloc[i][Scols.iloc[i]]
                        for i in np.arange(0, len(Scols))]
        readsIn['overlap'] = (readsIn['E'] - readsIn['S']
                              ) / readsIn['readLength']
        readsIn['sliceLength'] = (
            readsIn['overlap'] * readsIn['readLength']).astype(int)
        readsIn['outsideLength'] = readsIn['readLength'] - \
            readsIn['sliceLength']

        # select reads based on interval overlap filter
        if self.overlap == 0.0:
            readsIn = readsIn[readsIn['overlap'] > self.overlap]
        else:
            readsIn = readsIn[readsIn['overlap'] >= self.overlap]

        # choose methylation information that is only contained in the read overlap, or in the entire read
        if self.readSlice == False:
            s = 'meth'
        if self.readSlice == True:
            s = 'slice'
            # get slices of meth seq contained in overlap:
            slices = []
            for index, row in readsIn.iterrows():
                sliceStart = row['S'] - row['read_start']
                sliceEnd = len(row['meth']) - (row['read_end'] - row['E'])
                slices.append(row['meth'][sliceStart:sliceEnd])
            readsIn[s] = slices

        return readsIn, intervalname

    def stitchReads(self, overlapReads):
        
        if self.readSlice == False:
            s = 'meth'
        if self.readSlice == True:
            s = 'slice'

        stichedReads = []

        for g, d in overlapReads.groupby('read'):

            # if read pair
            if len(d) == 2:

                # if the read overlaps do not overlap themselves, stitch together the slices
                if (d['read_end'].iloc[0] < d['read_start'].iloc[1]):
                    stitchedMeth = d[s].iloc[0] + d[s].iloc[1]

                if (d['read_end'].iloc[1] < d['read_start'].iloc[0]):
                    stitchedMeth = d[s].iloc[1] + d[s].iloc[0]

                # if the read overlaps are identical:
                if (d['read_start'].iloc[0] == d['read_start'].iloc[1]) and (d['read_end'].iloc[0] == d['read_end'].iloc[1]):
                    stitchedMeth = d[s].iloc[0]

                # if end of R1 overlaps beginning of R2:
                if (d['read_end'].iloc[0] > d['read_start'].iloc[1]) and (d['read_start'].iloc[1] > d['read_start'].iloc[0]):
                    if s == 'slice':
                        readOverlap = d['E'].iloc[0] - d['S'].iloc[1]
                        stitchedMeth = d[s].iloc[0] + \
                            d[s].iloc[1][readOverlap:]
                    if s == 'meth':
                        readOverlap = d['read_end'].iloc[0] - \
                            d['read_start'].iloc[1]
                        stitchedMeth = d[s].iloc[0] + \
                            d[s].iloc[1][readOverlap:]

                # if end of R2 overlaps beginning of R1:
                if (d['read_end'].iloc[1] > d['read_start'].iloc[0]) and (d['read_start'].iloc[0] > d['read_start'].iloc[1]):
                    if s == 'slice':
                        readOverlap = d['E'].iloc[1] - d['S'].iloc[0]
                        stitchedMeth = d[s].iloc[1] + \
                            d[s].iloc[0][readOverlap:]
                    if s == 'meth':
                        readOverlap = d['read_end'].iloc[1] - \
                            d['read_start'].iloc[0]
                        stitchedMeth = d[s].iloc[1] + \
                            d[s].iloc[0][readOverlap:]

            # if single read
            if len(d) == 1:
                stitchedMeth = d[s].iloc[0]

            numU = stitchedMeth.count('z')
            numM = stitchedMeth.count('Z')

            # only keep reads that cover at least 1 CpG
            if numU + numM == 0:
                continue
            else:
                md = numM / float(numU + numM)
                chrom = d['chr'].iloc[0]
                if len(d) == 2:
                    if s == 'slice':
                        methSeqStart = np.min([d['S'].iloc[0], d['S'].iloc[1]])
                        methSeqStop = np.max([d['E'].iloc[0], d['E'].iloc[1]])
                    if s == 'meth':
                        methSeqStart = np.min(
                            [d['read_start'].iloc[0], d['read_start'].iloc[1]])
                        methSeqStop = np.max(
                            [d['read_end'].iloc[0], d['read_end'].iloc[1]])
                else:
                    if s == 'slice':
                        methSeqStart = d['S'].iloc[0]
                        methSeqStop = d['E'].iloc[0]
                    if s == 'meth':
                        methSeqStart = d['read_start'].iloc[0]
                        methSeqStop = d['read_end'].iloc[0]

                stichedReads.append([chrom, d['int_start'].iloc[0], d['int_stop'].iloc[0],
                                     methSeqStart, methSeqStop, len(stitchedMeth), numU, numM, md])

        return stichedReads

    def uniqueReadCounts(self, sampleReads):

        cols = ['chr', 'interval_start', 'interval_stop', 'methSeqStart',
                'methSeqStop', 'methSeqLength', 'numU', 'numM', 'MD']
        sample_reads_df = pd.DataFrame(sampleReads, columns=cols)
        sample_reads_df.set_index(cols, inplace=True)

        # get counts of unique reads
        unique_read_counts = dict(Counter(sample_reads_df.index.tolist()))

        # index becomes the read, and the column becomes the number of that particular read in the sample
        unique_read_counts_df = pd.DataFrame.from_dict(
            unique_read_counts, orient='index')

        return unique_read_counts_df, cols
    
    @property
    def densityTable(self):

        merged_df = pd.DataFrame()

        for f in self.files:
            print(' ')
            print('extracting reads from: ' + path_leaf(f))
            sample_name = path_leaf(f).split('.', 1)[0]

            sample_reads = []

            if find_executable('samtools') is not None:
                lines = self.samtoolsReadLines(f)
            # if samtools does not exist, and sample is sam:
            elif bool(re.search(r'sam$', f, re.IGNORECASE)):
                lines = self.manualReadLines(f)
            # if samtools does not exist, and sample is bam:
            elif bool(re.search(r'bam$', f, re.IGNORECASE)):
                print('`samtools view` needed to access reads in BAM file.')
                print('skipping file...')
                continue

            read_df = pd.DataFrame(
                lines, columns=['read', 'chr', 'read_start', 'meth'])
            # remove 'XM:Z:' at beginning of meth sequence
            read_df['meth'] = read_df['meth'].str.split(':', 2).str[2]
            read_df['readLength'] = read_df['meth'].str.len()
            read_df['read_start'] = read_df['read_start'].astype(int)
            read_df['read_end'] = read_df['read_start'] + read_df['readLength']

            for interval in self.intervals.itertuples():
                overlapReads, intervalName = self.overlappingReads(
                    read_df, interval)
                stitchedReads = self.stitchReads(overlapReads)
                if len(stitchedReads) == 0:
                    print('Warning: No reads found overlapping with ' +
                          intervalName + ' in ' + str(path_leaf(f)))
                    continue
                else:
                    print(' ' + str(len(stitchedReads)) +
                          ' unique sequences containing CpGs in ' + intervalName)
                    for i in stitchedReads:
                        sample_reads.append(i)

            # check to see if any reads were selected:
            if len(sample_reads) == 0:
                print(
                    'Warning: No reads found overlapping any intervals in: ' + str(path_leaf(f)))
                continue
            else:
                print(' ' + str(len(sample_reads)) +
                      ' total unique sequences containing CpGs')

                sampleUniqueReadCounts, cols = self.uniqueReadCounts(sample_reads)

                # if same sample name had been used, append a counter to the end
                if sample_name in merged_df.columns:
                    sample_name = sample_name + '.' + \
                        str(len(
                            [i for i in merged_df.columns.tolist() if sample_name in i]))

                sampleUniqueReadCounts.columns = [sample_name]

                # append the unique_read_counts_df sample dataframe to the final dataframe
                merged_df = merged_df.join(sampleUniqueReadCounts, how='outer')

        # cleanup merged_df
        if merged_df.empty == True:
            print(' No samples found in ' + str(self.fileDirectory) + '.')
            print(' Check to make sure path is correct or samples have proper extensions.')
            sys.exit(0)
        else:

            print(' ')
            print('Cleaning table...')
            # cleanup final dataframe
            sample_names = merged_df.columns.tolist()
            merged_df['chr'] = [i[0] for i in merged_df.index]
            merged_df['interval_start'] = [i[1] for i in merged_df.index]
            merged_df['interval_stop'] = [i[2] for i in merged_df.index]
            merged_df['methSeqStart'] = [i[3] for i in merged_df.index]
            merged_df['methSeqStop'] = [i[4] for i in merged_df.index]
            merged_df['methSeqLength'] = [i[5] for i in merged_df.index]
            merged_df['numU'] = [i[6] for i in merged_df.index]
            merged_df['numM'] = [i[7] for i in merged_df.index]
            merged_df['MD'] = [i[8] for i in merged_df.index]
            merged_df.reset_index(inplace=True)
            column_order = cols + sample_names
            merged_df = merged_df[column_order]
            merged_df.sort_values(
                ['chr', 'methSeqStart', 'methSeqStop'], inplace=True)

        return merged_df

# boxplot class
class boxplot(object):
    '''
    Each column is its own boxplot.
    '''

    @staticmethod
    def lighten_color(color, amount=0.3):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.

        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        import matplotlib.colors as mc
        import colorsys
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

    def __init__(self, df, colors=None, xtitle=None, ytitle=None, title=None):

        self.xtitle = xtitle
        self.ytitle = ytitle
        self.title = title

        if colors == None:
            colors = ['gray'] * len(df)
            colors = itertools.cycle(colors)
        else:
            colors = itertools.cycle(colors)

        labels = []
        data = []
        color = []

        # check if pandas dataframe, else convert to one
        if isinstance(df, pd.DataFrame):
            self.df = df
        elif isinstance(df, np.array):
            self.df = pd.DataFrame(df)
        elif isinstance(df, list):
            self.df = pd.DataFrame(df)

        self.boxRange = np.arange(len(self.df.columns))
        for l, v, c in zip(self.df.columns.tolist(), self.boxRange, [next(colors) for i in self.boxRange]):
            labels.append(l)
            data.append(self.df.iloc[:, v].dropna().tolist())
            color.append(c)

        self.labels = labels
        self.values = data
        self.colors = color
        self.size = (len(df), 4)
        self.positions = list(np.arange(1, len(self.labels)*2, 2))
        self.maxVal = max(
            [item for sublist in self.values for item in sublist])
        self.facecolors = [self.lighten_color(c) for c in self.colors]

    def plot(self):

        fig, ax = plt.subplots(figsize=self.size)
        bp = ax.boxplot(self.values, positions=self.positions,
                        widths=1, patch_artist=True)

        for flier in bp['fliers']:
            flier.set(marker='', color='black')
        for whisker in bp['whiskers']:
            whisker.set(color='black', linewidth=2)
        for cap in bp['caps']:
            cap.set(color='black', linewidth=2)
        for median in bp['medians']:
            median.set(color='black', linewidth=2)

        for i in self.boxRange:
            bp['boxes'][i].set(color=self.colors[i], linewidth=2, alpha=0.9)
            bp['boxes'][i].set(facecolor=self.facecolors[i])
            scatter = ax.scatter(x=np.random.normal(self.positions[i], 0.1, size=len(self.values[i])),
                                 y=self.values[i], c=self.colors[i], marker='.', edgecolors='', s=50, zorder=10)

        ax.set_xlim([0, max(self.positions)+1])
        ax.set_ylim([0, self.maxVal * 1.1])
        ax.set_ylabel(self.ytitle, fontsize=18)
        plt.yticks(np.linspace(0, self.maxVal*1.05, 5), ['0.0']+['%.1e' % i for i in np.linspace(0, self.maxVal*1.05, 5)[1:]],
                   fontsize=12)
        ax.set_xlabel(self.xtitle, fontsize=24)
        ax.set_xticklabels(self.labels, fontsize=12, rotation=45)
        plt.title(self.title, fontsize=18)
        return plt

# boxplot class specific for just case and controls comparison, with stats


class boxplot2sets(boxplot):
    '''
    Can select the case column and control column to compare as one boxplot in this subclass.
    If not indicated then case column is first column in df and control column is second.
    '''

    def __init__(self, df, colors=None, xtitle=None, ytitle=None, title=None,
                 case=None, control=None):

        # inherit args from parent boxplot class
        super(boxplot2sets, self).__init__(df, colors, xtitle, ytitle, title)

        self.labels = []
        # select cases/ctrl vals based on order in df (0,1) or columns labels if strings
        if case != None:
            self.case = df[case].dropna().tolist()
            self.labels.append(case)
        else:
            self.case = df.iloc[:, 0].dropna().tolist()
            self.labels.append('cases')
        if control != None:
            self.control = df[control].dropna().tolist()
            self.labels.append(control)
        else:
            self.control = df.iloc[:, 1].dropna().tolist()
            self.labels.append('controls')

        # change some attributes b/c comparing 2 sample sets only
        self.boxRange = np.arange(2)
        self.size = (2, 4)
        self.values = [self.case, self.control]
        self.positions = [1, 3]

    @property
    def ranksum(self):

        s, p = stats.ranksums(self.case, self.control)
        return p

    @property
    def cutoffVal(self):

        case_labels = [1 for i in self.case]
        ctrl_labels = [0 for i in self.control]
        roc_values = [item for sublist in self.values for item in sublist]
        roc_labels = case_labels + ctrl_labels

        fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
        roc_auc = auc(fpr, tpr)

        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        return optimal_threshold

    def plot(self, stats=False, thresh=False):

        plt = super(boxplot2sets, self).plot()

        if stats is not False:
            if 0.01 <= self.ranksum < 0.05:
                plt.text(x=1.7, y=self.maxVal * 1.10, s='*', fontsize=24)
            if 0.001 <= self.ranksum < 0.01:
                plt.text(x=1.6, y=self.maxVal * 1.10, s='**', fontsize=24)
            if self.ranksum < 0.001:
                plt.text(x=1.5, y=self.maxVal * 1.10, s='***', fontsize=24)
            if self.ranksum >= 0.05:
                plt.text(x=1.5, y=self.maxVal * 1.12, s='ns', fontsize=24)
            plt.title(self.title, fontsize=18, y=1.10)

        if thresh is not False:
            plt.axhline(y=self.cutoffVal, linestyle='--', color='k')

        return plt

# class for stacked barplots


class stackedBarplot():
    """
    df should have MD vals as indices and sample read counts as columns.
    For future/general use, values to colorby could be a column and not index
    """

    def __init__(self, df, colormap='coolwarm', colorby=None, columns=None,
                 xtitle=None, ytitle=None, title=None, colorbarlabel=None):

        self.xtitle = xtitle
        self.ytitle = ytitle
        self.title = title
        self.colorbarlabel = colorbarlabel

        # color data by values in a column or use the index
        if colorby != None:
            self.colorVals = sorted(list(set(df[colorby].values)))
        else:
            self.colorVals = sorted(list(set(df.index.values)))

        self.cmap = plt.get_cmap(colormap)
        self.colors = self.cmap(self.colorVals)

        # select sample columns to plot
        if columns != None:
            self.values = df[columns]
        else:
            self.values = df

        if isinstance(self.values, pd.DataFrame):
            self.nsamples = len(self.values.columns)
            self.maxval = max(list(self.values.sum()))
        if isinstance(self.values, pd.Series):
            self.nsamples = 1
            self.maxval = max(self.values)

    def plot(self):

        fig = plt.figure()

        # plot the data
        ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])
        plot = self.values.T.plot(kind='bar', stacked=True,
                                  ax=ax, color=self.colors, width=1, legend=False)
        plt.yticks(fontsize=14)
        plt.title(self.title, fontsize=24)

        ax.set_xlim([-0.5, self.nsamples - 0.5])
        ax.set_ylim(0, self.maxval)
        ax.set_ylabel(self.ytitle, fontsize=18)

        # make bottom x-axis labels the number of reads
        #ax.set_xticklabels([str(int(i)) for i in self.values.sum().values], rotation=45, fontsize=10)
        #ax.set_xlabel('reads covering locus', fontsize=14)

        # colorbar
        ax2 = fig.add_axes([0.1, 1.25, 0.5, 0.05])
        cb = mpl.colorbar.ColorbarBase(ax2, cmap=self.cmap, orientation='horizontal',
                                       norm=mpl.colors.Normalize(vmin=0, vmax=1))
        cb.set_label(self.colorbarlabel, fontsize=18)
        cb.ax.xaxis.set_ticks_position('top')
        #cb.set_clim(0.0, 1)
        # mpl.cm.ScalarMappable.set_clim([0.0,1.0])

        # x-axis label = sample type and number of samples
        ax4 = ax.twiny()
        ax4.set_xticks([])
        if isinstance(self.values, pd.DataFrame):
            ax4.set_xticklabels([self.values.columns], fontsize=12)
        if isinstance(self.values, pd.Series):
            ax4.set_xticklabels([self.values.name], fontsize=12)
        ax4.set_xlabel(self.xtitle, fontsize=24)
        ax4.xaxis.set_label_coords(0.5, -0.4)

        return plt

# class for histograms


class histogram():
    '''
    fractional distribution of reads in entire sample set with different MDs
    df indices should be MD values and columns are sample read counts.
    For future/general use, sample set can be all columns in df, or specific columns can be selected
    '''

    @staticmethod
    def lighten_color(color, amount=0.7):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.

        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        import matplotlib.colors as mc
        import colorsys
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

    def __init__(self, df, cases=None, controls=None, meCpGs=None,
                 casecolor='red', controlcolor='blue', defaultcolor='gray',
                 caselabel=None, controllabel=None, defaultlabel=None,
                 xtitle=None, ytitle=None, title=None, legend=False):

        self.cases = cases
        self.controls = controls

        self.casecolor = self.lighten_color(casecolor)
        self.controlcolor = self.lighten_color(controlcolor)
        self.defaultcolor = defaultcolor

        self.caselabel = caselabel
        self.controllabel = controllabel
        self.defaultlabel = defaultlabel

        self.xtitle = xtitle
        self.ytitle = ytitle
        self.title = title

        self.meCpGs = meCpGs

        self.legend = legend

        # sum of reads for each MD/index value depending on columns/samples sets used
        if isinstance(df, pd.DataFrame):

            if self.cases != None and self.controls != None:
                self.casevals = df[cases].sum(axis=1)
                self.controlvals = df[controls].sum(axis=1)
                self.values = list(self.casevals.values) + \
                    list(self.controlvals.values)

            if self.controls != None and self.cases == None:
                self.casevals = None
                self.controlvals = df[controls].sum(axis=1)
                self.values = list(self.controlvals.values)

            if self.controls == None and self.cases != None:
                self.casevals = df[cases].sum(axis=1)
                self.controlvals = None
                self.values = list(self.casevals.values)

            if self.cases == None and self.controls == None:
                self.casevals = None
                self.controlvals = None
                self.values = df.sum(axis=1)

        # series has one value per index
        if isinstance(df, pd.Series):

            self.casevals = None
            self.controlvals = None
            self.values = df

    @property
    def adjust(self):
        # in the case that the 'sums' for each index are actually fractions and not integers
        # then will need to convert to integers using this adjustment value.
        # could happen if sample read counts are adjusted by relative sample fractions.
        # make the lowest 'sum' an integer of at least 10 or greater.
        # using this may increase the number of occurances in list for hist
        # but also use this value to lower the occurance weights to equalize hist

        if isinstance(self.values, list):
            adjust = int(
                round(10.0/np.min([i for i in self.values if i != 0.0])+1))
            return adjust
        else:
            values = list(self.values.values)
            adjust = int(round(10.0/np.min([i for i in values if i != 0.0])+1))
            return adjust

    @property
    def caseDistr(self):

        # will be pd.Series if self.cases != None
        if isinstance(self.casevals, pd.Series):
            caseDistr = []
            for loc, val in enumerate(self.casevals):
                caseDistr += int(val*self.adjust) * \
                    [float(self.casevals.index[loc])]
        else:
            caseDistr = None

        return caseDistr

    @property
    def caseweights(self):

        # will be pd.Series if self.cases != None
        if isinstance(self.casevals, pd.Series):
            # weights = 1/num reads or occurances, ex: 100 reads --> 1/100 weight for each read
            caseweights = (np.ones_like(self.caseDistr) /
                           float(np.nansum(self.casevals)))/self.adjust
        else:
            caseweights = None

        return caseweights

    @property
    def controlDistr(self):

        # will be pd.Series if self.controls != None
        if isinstance(self.controlvals, pd.Series):
            controlDistr = []
            for loc, val in enumerate(self.controlvals):
                controlDistr += int(val*self.adjust) * \
                    [float(self.controlvals.index[loc])]
            controlDistr = controlDistr
        else:
            controlDistr = None

        return controlDistr

    @property
    def controlweights(self):

        # will be pd.Series if self.controls != None
        if isinstance(self.controlvals, pd.Series):
            # weights = 1/num reads or occurances, ex: 100 reads --> 1/100 weight for each read
            controlweights = (np.ones_like(self.controlDistr) /
                              float(np.nansum(self.controlvals)))/self.adjust
        else:
            controlweights = None

        return controlweights

    @property
    def valueDistr(self):

        # will be true if entire df to be used or a series is being used
        if isinstance(self.values, pd.Series) or (self.cases == None and self.controls == None):
            valueDistr = []
            for loc, val in enumerate(self.values):
                valueDistr += int(val*self.adjust) * \
                    [float(self.values.index[loc])]
        else:
            valueDistr = None

        return valueDistr

    @property
    def valueweights(self):

        # will be true if entire df to be used or a series is being used
        if isinstance(self.values, pd.Series) or (self.cases == None and self.controls == None):
            # weights = 1/num reads or occurances, ex: 100 reads --> 1/100 weight for each read
            valueweights = (np.ones_like(self.valueDistr) /
                            float(np.nansum(self.values)))/self.adjust
        else:
            valueweights = None

        return valueweights

    @property
    def bins(self):
        # bins (based on number of possible meCpGs in reads)
        # if # CpGs is less than 10, then would like each bin to represent reads of a single MD type.
        # for hist bin edges, by default: [1,2), [2,3), [3,4), [4,5]
        # so for 12.5% MD increments (8 CpGs):
        # [0, 0.1249), [0.1249, 0.2499), [0.2499, 0.3749), [0.3749, 0.4999), [0.4999, 0.6249), [0.6249, 0.7499), [0.7499, 0.8749), [0.8749, 1.0]
        # 0 MD reads,  0.125 MD reads,   0.25 MD reads,    0.375 MD reads,   0.50 MD reads,    0.675 MD reads,   0.75 MD reads,    0.875 and 1.0 MD reads
        # Otherwise, bins can just span segments of 10% MD increments.
        # will also default to 10% increments if no bins given

        if self.meCpGs == None:
            bins = 10
            return bins
        elif self.meCpGs >= 10:
            bins = 10
            return bins
        else:
            xtick_pos = list(np.linspace(
                0.0, 1.0, self.meCpGs, endpoint=False)) + [1.0]
            bins = [0.0] + [i-0.0001 for i in xtick_pos[1:-1]] + [1.0]
            return bins

    @property
    def xtick_pos(self):

        if self.meCpGs == None or self.meCpGs >= 10:
            xtick_pos = list(np.linspace(
                0.0, 1.0, self.bins, endpoint=False)) + [1.0]
        else:
            xtick_pos = list(np.linspace(
                0.0, 1.0, self.meCpGs, endpoint=False)) + [1.0]
        return xtick_pos

    def plot(self):

        f, ax = plt.subplots()

        if isinstance(self.caseDistr, list):
            cases_n, cases_bins, cases_patches = ax.hist(self.caseDistr,
                                                         bins=self.bins,
                                                         density=False,
                                                         weights=self.caseweights,
                                                         color=self.casecolor,
                                                         align='mid',
                                                         range=(0.0, 1.0),
                                                         alpha=0.7,
                                                         label=self.caselabel)
        if isinstance(self.controlDistr, list):
            ctrl_n, ctrl_bins, ctrl_patches = ax.hist(self.controlDistr,
                                                      bins=self.bins,
                                                      density=False,
                                                      weights=self.controlweights,
                                                      color=self.controlcolor,
                                                      align='mid',
                                                      range=(0.0, 1.0),
                                                      alpha=0.7,
                                                      label=self.controllabel)
        if isinstance(self.valueDistr, list):
            _n, _bins, _patches = ax.hist(self.valueDistr,
                                          bins=self.bins,
                                          density=False,
                                          weights=self.valueweights,
                                          color=self.defaultcolor,
                                          align='mid',
                                          range=(0.0, 1.0),
                                          alpha=0.7,
                                          label=self.defaultlabel)

        # focus on the methylated epialleles
        # background presumably mostly MD=0, so ignore and let go off axis
        # get list of bar heights except the <10% MD background bars and take max value as limit
        if isinstance(self.caseDistr, list) and isinstance(self.controlDistr, list):
            heights = list(cases_n[1:]) + list(ctrl_n[1:])

        if self.caseDistr == None and isinstance(self.controlDistr, list):
            heights = list(ctrl_n[1:])

        if isinstance(self.caseDistr, list) and self.controlDistr == None:
            heights = list(cases_n[1:])

        if isinstance(self.valueDistr, list):
            heights = list(_n[1:])

        plt.ylim([0, max(heights)*1.25])
        ytick_pos = np.linspace(0.0, max(heights)*1.20, 5)

        # set bins on x-axis:
        plt.xticks(self.xtick_pos, ['0%'] + [str(round(md*100, 1))+'%' for md in self.xtick_pos[1:]],
                   fontsize=18, rotation=30)

        # set range on y-axis:
        plt.yticks(ytick_pos, ['0.0'] + ['%.1e' %
                                         i for i in ytick_pos[1:]], fontsize=18)

        ax.set_title(self.title, fontsize=22)
        plt.ylabel(self.ytitle, fontsize=20, labelpad=40)
        plt.xlabel(self.xtitle, fontsize=20, rotation=0)
        ax.yaxis.set_label_coords(-0.23, 0.55)

        if self.legend == True:
            plt.legend(loc="upper center", fontsize=14,
                       edgecolor='k', bbox_to_anchor=(0.55, 0.90))

        return plt

# class for heatmaps


class heatmap():

    def __init__(self, matrix, xticks=None, yticks=None, colormap='hot',
                 xtitle=None, ytitle=None, title=None):

        self.matrix = matrix

        self.xtitle = xtitle
        self.ytitle = ytitle
        self.title = title

        self.xtickpos = np.arange(0, self.matrix.shape[1], 10)
        self.ytickpos = np.arange(0, self.matrix.shape[0], 4)

        if isinstance(xticks, np.ndarray):
            self.xticks = xticks
        elif isinstance(xticks, list):
            self.xticks = np.array(xticks)
        elif xticks == None:
            self.xticks = np.arange(0, self.matrix.shape[1])
        else:
            print('xticks needs to be 1D numpy.ndarray or list type')
            self.xticks = np.arange(0, self.matrix.shape[1])

        if isinstance(xticks, np.ndarray):
            self.yticks = list(yticks)
        elif isinstance(yticks, list):
            self.yticks = yticks
        elif yticks == None:
            self.yticks = list(np.arange(0, self.matrix.shape[0]))
        else:
            print('yticks needs to be 1D numpy.ndarray or list type')
            self.yticks = list(np.arange(0, self.matrix.shape[0]))

        # if yticks == None:
        #     self.yticks = np.arange(0, self.matrix.shape[0])
        # else:
        #     self.yticks = yticks

        self.cmap = plt.get_cmap(colormap)
        #self.colors = self.cmap(self.colorVals)

    def plot(self):

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])
        ax.imshow(self.matrix, cmap=self.cmap, vmin=0,
                  vmax=1)  # option: remove vmin, vmax

        plt.xticks(self.xtickpos, ['0.0'] + ['%.1e' % x for x in self.xticks[self.xtickpos][1:]],
                   rotation=45, fontsize=14)
        plt.yticks(self.ytickpos, [str(y) for y in np.array(self.yticks[::-1])[self.ytickpos]],
                   fontsize=16)

        plt.xlabel(self.xtitle, fontsize=20)
        plt.ylabel(self.ytitle, fontsize=20, labelpad=50, rotation=0)
        ax.yaxis.set_label_coords(-0.26, 0.3)

        # colorbar
        ax2 = fig.add_axes([0.30, 0.77, 0.5, 0.03])
        cb = mpl.colorbar.ColorbarBase(ax2, cmap=self.cmap, orientation='horizontal',
                                       norm=mpl.colors.Normalize(vmin=0, vmax=1))  # option: vmax = np.amax(tpr_matrix)
        cb.ax.set_title(self.title, fontsize=24)
        cb.ax.tick_params(labelsize=14)
        cb.ax.xaxis.set_ticks_position('bottom')
        # cb.set_clim(0.0, 1.0) # option: 0, np.amax(tpr_matrix)

        return plt

# class for ROC curves


class rocplot():
    '''
    Dataframe should have two columns, one contaning values for cases and other for controls.
    '''

    def __init__(self, df, cases=None, controls=None):

        import pandas as pd
        import numpy as np
        from sklearn.metrics import roc_curve

        # select cases/ctrl vals based on order in df (0,1) or columns labels if strings
        if cases != None:
            self.caseVals = df[cases].dropna().tolist()
        else:
            self.caseVals = df.iloc[:, 0].dropna().tolist()
        if controls != None:
            self.controlVals = df[controls].dropna().tolist()
        else:
            self.controlVals = df.iloc[:, 1].dropna().tolist()

        self.values = self.caseVals + self.controlVals
        self.labels = [1 for i in self.caseVals] + \
            [0 for i in self.controlVals]
        self.fpr, self.tpr, self.thresholds = roc_curve(
            self.labels, self.values)

        optimal_idx = np.argmax(self.tpr - self.fpr)
        self.optimal_threshold = self.thresholds[optimal_idx]
        self.sensitivity = self.tpr[optimal_idx]
        self.specificity = 1 - self.fpr[optimal_idx]

    @property
    def AUC(self):

        roc_auc = auc(self.fpr, self.tpr)
        return roc_auc

    def plot(self):

        lw = 2

        plt.plot(self.fpr, self.tpr, color='red', lw=lw,
                 label='MDBC AUC = %0.2f' % self.AUC)

        plt.plot([0, 1], [0, 1], color='k', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.00])
        plt.ylim([0.0, 1.00])
        plt.xlabel('FPR', fontsize=26)
        plt.ylabel('TPR', fontsize=26)
        plt.xticks(np.arange(0, 1.1, .2), [str(round(i, 2))
                                           for i in np.arange(0, 1.1, .2)], fontsize=20)
        plt.yticks(np.arange(0, 1.1, .2), [str(round(i, 2))
                                           for i in np.arange(0, 1.1, .2)], fontsize=20)
        plt.legend(loc="lower right", fontsize=12, edgecolor='k')

        return plt

# class for manipulating methylation density tables for input into plotting scripts


class mdbc():

    def __init__(self, df, cases, controls, fractions=None, mdcutoffs=None, hdf_label='mdbc.mdvals.h5'):

        self.densityTable = pd.read_csv(df)
        self.mdvalues = sorted(list(set(self.densityTable['MD'].values)))

        if isinstance(mdcutoffs, np.ndarray):
            self.mdcutoffs = mdcutoffs
        elif isinstance(mdcutoffs, list):
            self.mdcutoffs = mdcutoffs
        elif mdcutoffs == None:
            # methylation density (md) cutoffs (0%, 5%, 10%, etc...) to 100%
            self.mdcutoffs = [round(i,2) for i in np.arange(0.0, 1.05, 0.05)]

        if len(cases) == 1:
            ID = cases[0]
            self.cases = [col for col in self.densityTable.columns if str(ID) in col]
        else:
            self.cases = cases
        if len(controls) == 1:
            ID = controls[0]
            self.controls = [col for col in self.densityTable.columns if str(ID) in col]
        else:
            self.controls = controls

        self.samples = list(self.cases) + list(self.controls)

        self.sampleMethylatedReadCounts = self.densityTable[['MD'] + self.samples].set_index(
            'MD').loc[[i for i in set(self.mdvalues) if i != 0]].sum()
        self.sampleTotalReadCounts = self.densityTable[self.samples].sum()
        self.noValueSamples = list(
            self.sampleTotalReadCounts[self.sampleTotalReadCounts == 0].index)

        self.cases_filt = [i for i in self.cases if i not in self.noValueSamples]
        self.controls_filt = [
            i for i in self.controls if i not in self.noValueSamples]
        self.casesAndControls_filt = self.cases_filt + self.controls_filt
        self.sampleMethylatedReadCounts_filt = self.sampleMethylatedReadCounts[~self.sampleMethylatedReadCounts.index.isin(
            self.noValueSamples)]
        self.sampleTotalReadCounts_filt = self.sampleTotalReadCounts[~self.sampleTotalReadCounts.index.isin(
            self.noValueSamples)]

        if isinstance(fractions, pd.DataFrame):
            self.inputFracs = dict(
                zip([str(i) for i in fractions['samples'].values], fractions['fractions'].values))
            self.sampleFracs = [float(self.inputFracs[i])
                                for i in self.casesAndControls_filt]
            self.casesFracs = [float(self.inputFracs[i])
                               for i in self.cases_filt]
            self.ctrlFracs = [float(self.inputFracs[i])
                              for i in self.controls_filt]
        else:
            self.inputFracs = 1.0
            self.sampleFracs = [1.0 for i in self.casesAndControls_filt]
            self.casesFracs = [1.0 for i in self.cases_filt]
            self.ctrlFracs = [1.0 for i in self.controls_filt]
        
        self.hdf_label = hdf_label

    @property
    def sampleMethReadCountsAdj(self):
        '''
        Returns normalized table (based on relative sample fractions)
        of the sample count of all methylated reads.
        First column for cases and second column for controls.
        '''

        # note that for MD cutoff of 0, reads of MD > 0 used, so only methylated reads are of interest.
        # Therefore, just using total reads (which would include reads with MD = 0) not needed to set thresh
        # and also samples not classified using just total amount of reads recorded

        sampleMethReadCountsAdj = self.sampleMethylatedReadCounts_filt / self.sampleFracs
        cases = sampleMethReadCountsAdj.loc[self.cases_filt].values.tolist()
        controls = sampleMethReadCountsAdj.loc[self.controls_filt].values.tolist(
        )
        df = pd.DataFrame([cases, controls]).T
        df.columns = ['cases', 'controls']

        return df

    @property
    def sampleMethReadEFsAdj(self):
        '''
        Returns normalized table (based on relative sample fractions)
        of the sample fraction of all methylated reads.
        '''

        sampleMethReadEFsAdj = (self.sampleMethylatedReadCounts_filt /
                                self.sampleTotalReadCounts_filt) / self.sampleFracs
        cases = sampleMethReadEFsAdj.loc[self.cases_filt].values.tolist()
        controls = sampleMethReadEFsAdj.loc[self.controls_filt].values.tolist()
        df = pd.DataFrame([cases, controls]).T
        df.columns = ['cases', 'controls']

        return df

    @property
    def CpGcolumns(self):
        '''
        Table with meC, unmeC and total C counts for each read (row) in densityTable.
        '''

        MD_cols = self.densityTable[['numU', 'numM', 'MD']]
        MD_cols2 = MD_cols.copy()
        MD_cols2['C'] = MD_cols[['numU', 'numM']].sum(axis=1)

        return MD_cols2

    @property
    def normalizedDensTable(self):
        '''
        Normalize the density table to the fraction of each sample analyzed.
        Fractions could be relative amounst of each sample or also amount(mg, ex)/volume (uL, ex) loaded per sample
        '''

        return self.densityTable[self.casesAndControls_filt].div(self.sampleFracs)

    @property
    def sampleAveMeth(self):
        '''
        Returns table of average methylation values for each sample.
        '''

        # list to collect sample average methylation values:
        AveMeth = []

        df = pd.concat([self.CpGcolumns, self.normalizedDensTable], axis=1)

        for sample in self.casesAndControls_filt:
            # temporary table of meC and total C counts for each read
            temp = self.CpGcolumns.copy()
            counts = df[sample]  # read counts for each MD of the given sample
            # counts of C's covered by each read in sample
            temp['sample_C'] = temp['C'].values * counts.values
            total_C = np.nansum(temp['sample_C'])

            # compute average methylation (use to sort samples in plot)
            # (total number of meC's covered by reads / total C's covered by reads)
            ave_meC = np.nansum(temp['numM'].values *
                                counts.values) / float(total_C)
            AveMeth.append(ave_meC)

        # reorder samples based on average methylation values
        reorder = zip(AveMeth, self.casesAndControls_filt)
        reorder = sorted(reorder, key=lambda x: x[0])

        return pd.DataFrame(reorder, columns=['average methylation', 'samples'])

    @property
    def readCountsPerMDtables(self):
        '''
        Returns two dataframes, first is for filtered cases and second is for filtered controls.

        Dataframe contains the count of reads with a given MD (rows) for each sample (columns).
        '''

        # here, don't care about counts of unmethylated reads
        #sorted_MDs = [ i for i in self.mdvalues if i != 0.0 ]
        sorted_MDs = [i for i in self.mdvalues]

        # make new dataframe to populate with relative counts of methylated reads for each sample:
        counts_per_MD = pd.DataFrame(index=sorted_MDs)

        df = pd.concat([self.CpGcolumns, self.normalizedDensTable], axis=1)

        # sample order based on average methylation values
        samples = self.sampleAveMeth['samples'].values.tolist()

        for sample in samples:
            counts = []  # populate list with counts of reads for each methylation density
            for m in sorted_MDs:
                # get counts of reads with given MD
                selection = df[df['MD'] == m][sample].values
                counts.append(sum(selection))
            # append list of read counts for each MD for given sample to dataframe
            counts_per_MD[sample] = counts

        cases_counts_per_MD = counts_per_MD[[
            i for i in samples if i in self.cases_filt]]
        control_counts_per_MD = counts_per_MD[[
            i for i in samples if i in self.controls_filt]]

        return cases_counts_per_MD, control_counts_per_MD

    @property
    def readEFsPerMDtables(self):
        '''
        Returns two dataframes, first is for filtered cases and second is for filtered controls.

        Dataframe contains the fraction of reads with a given MD (rows) for each sample (columns).
        '''

        # make new dataframe to populate with relative fractions of reads of each MD type for each sample:
        # here, fractions of reads for each MD should add up to 1.
        # contribution of each read to given fraction is weighted by number of C's covered by the read
        fractions_per_MD = pd.DataFrame(index=self.mdvalues)

        df = pd.concat([self.CpGcolumns, self.normalizedDensTable], axis=1)

        # sample order based on average methylation values
        samples = self.sampleAveMeth['samples'].values.tolist()

        for sample in samples:
            # temporary table of meC and total C information of read patterns
            temp = self.CpGcolumns.copy()
            counts = df[sample]  # read counts for each MD of the given sample
            # counts of C's covered by each read in sample
            temp['sample_C'] = temp['C'].values * counts.values
            # total C's covered by all reads in sample
            total_C = np.nansum(temp['sample_C'])

            # weighted fraction
            # (number of C's covered by reads with given methylation density / total C's covered by reads in case)
            weighted_fracs = []
            for m in self.mdvalues:
                # get rows (read patterns) that have the given methylation density
                selection = temp[temp['MD'] == m]
                # weighted fraction of reads with given MD
                frac = np.nansum(selection['sample_C']) / float(total_C)
                weighted_fracs.append(frac)

            fractions_per_MD[sample] = weighted_fracs

        cases_fractions_per_MD = fractions_per_MD[[
            i for i in samples if i in self.cases_filt]]
        control_fractions_per_MD = fractions_per_MD[[
            i for i in samples if i in self.controls_filt]]

        return cases_fractions_per_MD, control_fractions_per_MD

    def sampleValsForMD(self, mdcutoff):
        '''
        Get sample read values for a given methylation density (MD) cutoff.

        Values are based on either the normalized counts of reads at or above the MD cutoff,
        or the relative sample fraction of reads at or above the MD cutoff.

        mdcutoff should be a float that is <= 1.0 and >= 0.0.
        '''

        if mdcutoff == 0.0:
            mdVals = [i for i in self.mdvalues if i > 0.0]
        else:
            mdVals = [i for i in self.mdvalues if i >= mdcutoff]

        countCasesTable = self.readCountsPerMDtables[0]
        countControlsTable = self.readCountsPerMDtables[1]

        efCasesTable = self.readEFsPerMDtables[0]
        efControlsTable = self.readEFsPerMDtables[1]

        caseCounts = []
        caseEFs = []
        for case in self.cases_filt:
            caseCounts.append(countCasesTable.loc[mdVals, case].sum())
            caseEFs.append(efCasesTable.loc[mdVals, case].sum())

        ctrlCounts = []
        ctrlEFs = []
        for ctrl in self.controls_filt:
            ctrlCounts.append(countControlsTable.loc[mdVals, ctrl].sum())
            ctrlEFs.append(efControlsTable.loc[mdVals, ctrl].sum())

        countVals = pd.DataFrame(
            {'cases': pd.Series(caseCounts), 'controls': pd.Series(ctrlCounts)})
        EFVals = pd.DataFrame(
            {'cases': pd.Series(caseEFs), 'controls': pd.Series(ctrlEFs)})

        return countVals, EFVals
    
    @property
    def storeSampleValuesPerMD(self):
        
        MDvalueKeys = {}
        for m in self.mdcutoffs:
            
            mdlabel = "_" + "_".join(str(m).split('.'))
            keys = [mdlabel + '_counts', mdlabel + '_efs']
            MDvalueKeys[m] = keys
            
            with pd.HDFStore(self.hdf_label) as store:
                if keys[0] in store:
                    pass
                else:
                    print(' Appending sample read count and fraction values for MD cutoff {} to HDF5 file {}'.format(m, path_leaf(self.hdf_label)))
                    vals = self.sampleValsForMD(m)
                    store[keys[0]] = vals[0]
                    store[keys[1]] = vals[1]
                    
        return MDvalueKeys

    @property
    def readCountCutoffRange(self):
        '''
        Define the range of possible read cutoffs to use for the number of methylated reads in a sample
        to classiy it as positive, given an MD cutoff.
        '''

        max_ef = pd.concat(self.readCountsPerMDtables,
                           axis=1).max().max() * 1.10

        return np.array(list(np.linspace(0.000, max_ef, 100, endpoint=False)) + [max_ef])

    @property
    def readEFCutoffRange(self):
        '''
        Define the range of possible read cutoffs to use for the fraction of methylated reads in a sample
        to classiy it as positive, given a MD cutoff.
        '''

        max_ef = pd.concat(self.readEFsPerMDtables, axis=1).max().max() * 1.10
        # EFs will be fractions, limit to 1.0
        if max_ef > 1.0:
            return np.array(list(np.linspace(0.000, 1.0, 100, endpoint=False)) + [1.0])
        else:
            return np.array(list(np.linspace(0.000, max_ef, 100, endpoint=False)) + [max_ef])

    def buildMatrices(self, read_metric='count'):
        '''
        Generate TPR/FPR/TPR-FPR matrices for each MD and
        sample read count or fraction cutoff combination.

        MD cutoffs are methylation density cutoffs.

        Read cutoffs are the sample read counts or fractions at or above the given MD cutoff
        to call sample positive.
        '''

        if read_metric == 'count':
            i = 0
            cutoffRange = self.readCountCutoffRange
        if read_metric == 'fraction':
            i = 1
            cutoffRange = self.readEFCutoffRange

        tpr_matrix = []
        fpr_matrix = []

        for m in self.mdcutoffs:
            
            # this class property is a dict where MD vals are keys.
            # dict values are list of keys to HDF5 table of sample values
            # first key in list is for count values, second is for fraction values
            # this property also generates the HDF5 entries if they have not been produced yet
            key = self.storeSampleValuesPerMD[m][i]
            df = pd.read_hdf(self.hdf_label, key=key)
            
            case_vals = df['cases'].dropna().tolist()
            case_labels = [1 for i in case_vals]
            ctrl_vals = df['controls'].dropna().tolist()
            ctrl_labels = [0 for i in ctrl_vals]
            roc_values = case_vals + ctrl_vals
            roc_labels = case_labels + ctrl_labels
            roc_df = pd.DataFrame({'label': roc_labels, 'values': roc_values})

            tprs = []
            fprs = []

            for cutoff in cutoffRange:
                true_positives = len(
                    roc_df[(roc_df['label'] == 1) & (roc_df['values'] > cutoff)].index)
                false_negatives = len(
                    roc_df[(roc_df['label'] == 1) & (roc_df['values'] <= cutoff)].index)
                true_negatives = len(
                    roc_df[(roc_df['label'] == 0) & (roc_df['values'] <= cutoff)].index)
                false_positives = len(
                    roc_df[(roc_df['label'] == 0) & (roc_df['values'] > cutoff)].index)
                tpr = float(true_positives) / \
                    (float(true_positives) + float(false_negatives))
                fpr = 1.0 - (float(true_negatives) /
                             (float(true_negatives) + float(false_positives)))
                tprs.append(tpr)
                fprs.append(fpr)

            tpr_matrix.append(tprs)
            fpr_matrix.append(fprs)

        tpr_matrix = np.flipud(np.array(tpr_matrix))
        fpr_matrix = np.flipud(np.array(fpr_matrix))
        diff_matrix = tpr_matrix - fpr_matrix

        return tpr_matrix, fpr_matrix, diff_matrix
    
    def cutoffPerformance(self, read_metric='count'):
        '''
        Generate summary table showing for each MD cutoff, the optimal read cutoff
        and corresponding TPR, FPR, and AUC values.
        '''

        if read_metric == 'count':
            label = 'optimal read count cutoff'
            i = 0
        if read_metric == 'fraction':
            label = 'optimal read fraction cutoff'
            i = 1

        summary = pd.DataFrame()
        aucs = []
        optimalMetricCutoffs = []
        TPRs = []
        specificity = []

        for m in self.mdcutoffs:
            
            # this class property is a dict where MD vals are keys.
            # dict values are list of keys to HDF5 table of sample values
            # first key in list is for count values, second is for fraction values
            # this property also generates the HDF5 entries if they have not been produced yet
            key = self.storeSampleValuesPerMD[m][i]
            df = pd.read_hdf(self.hdf_label, key=key)
            
            case_vals = df['cases'].dropna().tolist()
            case_labels = [1 for i in case_vals]
            ctrl_vals = df['controls'].dropna().tolist()
            ctrl_labels = [0 for i in ctrl_vals]
            roc_values = case_vals + ctrl_vals
            roc_labels = case_labels + ctrl_labels
            roc_df = pd.DataFrame({'label': roc_labels, 'values': roc_values})

            fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
            roc_auc = auc(fpr, tpr)
            optimal_idx = np.argmax(tpr - fpr)

            aucs.append(roc_auc)
            optimalMetricCutoffs.append(thresholds[optimal_idx])
            TPRs.append(tpr[optimal_idx])
            specificity.append(1 - fpr[optimal_idx])

        summary['MD Cutoffs'] = self.mdcutoffs
        summary[label] = optimalMetricCutoffs
        summary['AUC'] = aucs
        summary['TPR'] = TPRs
        summary['1 - FPR'] = specificity

        return summary

    def optimalMDcutoff(self, summary):
        '''
        Return optimal MD for all MD cutoffs.
        Defined as the cutoff that has the largest positive difference of TPR - FPR
        for a particular read cutoff value.
        '''

        # optimal MD defined as the largest positive difference of TPR - FPR for all cutoffs tested:
        diff_vals = summary['TPR'].values - (1 - summary['1 - FPR'].values)
        max_diff = np.max(diff_vals)
        # in the case of ties, pick MD that also has highest AUC. If still a tie, then pick lowest MD cutoff.
        max_idx = summary.iloc[np.argwhere(
            diff_vals == max_diff).flatten()]['AUC'].idxmax()
        opt_md = summary.iloc[max_idx]['MD Cutoffs']

        return opt_md

    # def __str__():

    # def __repr__():


# Make way to also reverse for hypomethylation? So instead of looping through >=MD, its <= MD?

# main script:


def main(args=None):

    args, args_dict = get_arguments(args)

    if args_dict['cmd'] == None:
        get_arguments(' ')

    print('#--------------------------------------------------------------------------')
    print('#--------------------------------------------------------------------------')
    print('Time:' + str(timestamp))
    print('#--------------------------------------------------------------------------')
    print('#--------------------------------------------------------------------------')
    print('Command:')
    print(sys.argv)
    print('Args:')
    print(args)
    print('#--------------------------------------------------------------------------')

    if args_dict['cmd'] == 'DREAMtoMD':
        rawDreamingData = args_dict['input']
        tempsToMDs = args_dict['tempsToMDs']
        numberCpGs = args_dict['numberCpGs']
        meltResolution = args_dict['tempResolution']
        includeBackground = args_dict['input2bg']
        poissonAdjustment = args_dict['poisson']
        includedRows = args_dict['includedRows']
        outdir = args_dict['output']
        inputName = path_leaf(rawDreamingData)

        if outdir != None:
            if outdir.endswith('/'):
                pass
            else:
                outdir = outdir + '/'
            filename = outdir + 'DREAMtoMD.DT.' + inputName + '.' + timestamp + '.csv'
        else:
            filename = cwd + '/DREAMtoMD.DT.' + inputName + '.' + timestamp + '.csv'

        print('DREAMtoMD:')
        print('Converting {} to methylation density table {}.'.format(inputName, path_leaf(filename)))

        dream = dreamingToDensityTable(rawDreamingData, tempsToMDs, numberCpGs,
                                       meltResolution=meltResolution, includeBackground=includeBackground,
                                       poissonAdjustment=poissonAdjustment, includedRows=includedRows)

        densityTable = dream.densityTable()
        densityTable.to_csv(filename, index=False)

        print(' ')
        print('{} created.'.format(filename))
        print(' ')
        print('DREAMtoMD complete.')
        print(' ')

    if args_dict['cmd'] == 'READtoMD':
        fileDirectory = args_dict['input']
        intervals = args_dict['intervals']
        fileType = args_dict['fileType']
        overlap = args_dict['overlap']
        readSlice = args_dict['slice']
        outdir = args_dict['output']

        if outdir != None:
            if outdir.endswith('/'):
                pass
            else:
                outdir = outdir + '/'
            filename = outdir + 'READtoMD.DT.' + timestamp + '.csv'
        else:
            filename = cwd + '/READtoMD.DT.' + timestamp + '.csv'

        print('READtoMD:')
        print('Converting files in {} to methylation density table {}.'.format(fileDirectory, path_leaf(filename)))

        reads = readsToDensityTable(fileDirectory=fileDirectory, intervals=intervals,
                                    fileType=fileType, readSlice=readSlice, overlap=overlap)

        densityTable = reads.densityTable
        densityTable.to_csv(filename, index=False)

        print(' ')
        print('{} created.'.format(filename))
        print(' ')
        print('READtoMD complete.')
        print(' ')

    if args_dict['cmd'] == 'MDBC':
        df = args_dict['input']
        cases = args_dict['cases']
        controls = args_dict['controls']
        fractions = args_dict['fractions']
        mdcutoffs = args_dict['MDcutoffs']
        outdir = args_dict['output']
        hdf_label = args_dict['hdf_label']

        dfName = path_leaf(df).split('.csv')[0]

        if outdir != None:
            if outdir.endswith('/'):
                pass
            else:
                outdir = outdir + '/'
            filename = outdir + 'MDBC.' + dfName
        else:
            filename = cwd + '/MDBC.' + dfName

        if hdf_label != None:
            pass
        else:
            hdf_label = cwd + '/MDBC_H5DF.' + dfName + '.h5'

        classifier = mdbc(df=df, cases=cases, controls=controls,
                         fractions=fractions, mdcutoffs=mdcutoffs,
                         hdf_label=hdf_label)

        # Default outputs:
        print('MDBC Analysis')
        print('')

        if args_dict['ignoreCountsummary'] == False:
            print(' Generating performance summary for each MD using normalized read counts...')
            countPerformance = classifier.cutoffPerformance(read_metric='count')
            countPerformance.to_csv(filename + '.READ-COUNT-SUMMARY.csv')
            countOptMD = classifier.optimalMDcutoff(countPerformance)
            print(' Optimal MD cutoff (read counts) = {}'.format(countOptMD))
            countOptMDVals = classifier.sampleValsForMD(countOptMD)[0]
            countBox = boxplot2sets(df=countOptMDVals, colors=['red', 'blue'],
                ytitle='normalized read count', title='MDC = {}'.format(countOptMD))
            print(' p-val (cases vs controls) = {}'.format(countBox.ranksum))
            countROC = rocplot(countOptMDVals)
            print(' AUC = {}'.format(countROC.AUC))
            countBox.plot(stats=True).savefig(
                filename + '.READ-COUNT-OPTMD-BOX.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            countROC.plot().savefig(
                filename + '.READ-COUNT-OPTMD-ROC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            print('')

        if args_dict['ignoreEFsummary'] == False:
            print(' Generating performance summary for each MD using normalized read fractions...')
            efPerformance = classifier.cutoffPerformance(read_metric='fraction')
            efPerformance.to_csv(filename + '.EF-SUMMARY.csv')
            efOptMD = classifier.optimalMDcutoff(efPerformance)
            print(' Optimal MD cutoff (read fractions) = {}'.format(efOptMD))
            efOptMDVals = classifier.sampleValsForMD(efOptMD)[1]
            efBox = boxplot2sets(df=efOptMDVals, colors=['red', 'blue'],
                ytitle='normalized sample read fraction', title='MDC = {}'.format(efOptMD))
            print(' p-val (cases vs controls) = {}'.format(efBox.ranksum))
            efROC = rocplot(efOptMDVals)
            print(' AUC = {}'.format(efROC.AUC))
            efBox.plot(stats=True).savefig(
                filename + '.EF-OPTMD-BOX.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            efROC.plot().savefig(
                filename + '.EF-OPTMD-ROC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            print('')

        if args_dict['totalreadcounts'] == True:
            print(' Returning total methylated read counts csv: ' + path_leaf(filename) + '.METH-TOTREADS.csv')
            methCounts = classifier.sampleMethReadCountsAdj
            methCounts.to_csv(filename + '.METH-TOTREADS.csv')

        if args_dict['totalreadcountsPlot'] == True:
            print(' Returning total methylated read counts boxplot: ' + path_leaf(filename) + '.METH-TOTREADS.png')
            methReadCountsPlot = boxplot2sets(
                df=classifier.sampleMethReadCountsAdj, colors=['red', 'blue'],
                ytitle='normalized read count', title='methylated reads')
            print(' p-val (cases vs controls) = {}'.format(methReadCountsPlot.ranksum))
            methReadCountsPlot.plot(stats=True).savefig(
                filename + '.METH-TOTREADS.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['totalEf'] == True:
            print(' Returning total methylated read fraction csv: ' + path_leaf(filename) + '.METH-TOTEFs.csv')
            methEFs = classifier.sampleMethReadEFsAdj
            methEFs.to_csv(filename + '.METH-TOTEFs.csv')

        if args_dict['totalEfPlot'] == True:
            print(' Returning total methylated read fraction boxplot: ' + path_leaf(filename) + '.METH-TOTEFs.png')
            methReadEFsPlot = boxplot2sets(
                df=classifier.sampleMethReadEFsAdj, colors=['red', 'blue'],
                ytitle='normalized sample read fraction', title='methylated reads')
            print(' p-val (cases vs controls) = {}'.format(methReadEFsPlot.ranksum))
            methReadEFsPlot.plot(stats=True).savefig(
                filename + '.METH-TOTEFs.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['readcountsEachMD'] == True:
            print(' Returning read counts for each MD cutoff csv, cases and controls: ' + path_leaf(filename) + '.*-READS.csv')
            casesCounts, controlCounts = classifier.readCountsPerMDtables
            casesCounts.to_csv(filename + '.CASES-READS.csv')
            controlCounts.to_csv(filename + '.CONTROLS-READS.csv')

        if args_dict['casesReadCountPlot'] == True:
            print(' Returning cases read counts for each MD barplot: ' + path_leaf(filename) + '.CASES-READS-COUNTS.png')
            casesCounts, controlCounts = classifier.readCountsPerMDtables
            casesCountsBar = stackedBarplot(casesCounts, ytitle='normalized read counts',
                xtitle='cases', colorbarlabel='read methylation density')
            casesCountsBar.plot().savefig(
                filename + '.CASES-READ-COUNTS.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['controlsReadCountPlot'] == True:
            print(' Returning controls read counts for each MD barplot: ' + path_leaf(filename) + '.CONTROLS-READS-COUNTS.png')
            casesCounts, controlCounts = classifier.readCountsPerMDtables
            controlCountsBar = stackedBarplot(controlCounts, ytitle='normalized read counts',
                xtitle='controls', colorbarlabel='read methylation density')
            controlCountsBar.plot().savefig(
                filename + '.CONTROLS-READ-COUNTS.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['EfEachMD'] == True:
            print(' Returning read fractions for each MD csv, cases and controls: ' + path_leaf(filename) + '.*-EFS.png')
            casesEFs, controlEFs = classifier.readEFsPerMDtables
            casesEFs.to_csv(filename + '.CASES-EFS.csv')
            controlEFs.to_csv(filename + '.CONTROLS-EFS.csv')

        if args_dict['casesEfPlot'] == True:
            print(' Returning cases read fractions for each MD barplot: ' + path_leaf(filename) + '.CASES-READ-FRACs.png')
            casesEFs, controlEFs = classifier.readEFsPerMDtables
            casesEFBar = stackedBarplot(casesEFs, ytitle='normalized sample read fraction',
                xtitle='cases', colorbarlabel='read methylation density')
            casesEFBar.plot().savefig(
                filename + '.CASES-READ-FRACs.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['controlsEfPlot'] == True:
            print(' Returning controls read fractions for each MD barplot: ' + path_leaf(filename) + '.CONTROLS-READ-FRACs.png')
            casesEFs, controlEFs = classifier.readEFsPerMDtables
            controlEFBar = stackedBarplot(controlEFs, ytitle='normalized sample read fraction',
                xtitle='controls', colorbarlabel='read methylation density')
            controlEFBar.plot().savefig(
                filename + '.CONTROLS-READ-FRACs.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['readcountDistributionPlot'] == True:
            print(' Returning distribution of reads for each MD histogram: ' + path_leaf(filename) + '.READ-DISTR.png')
            combinedCounts = pd.concat(classifier.readCountsPerMDtables, axis=1)
            countsHist = histogram(combinedCounts, cases=classifier.cases_filt,
                                   controls=classifier.controls_filt, meCpGs=len(classifier.CpGcolumns),
                                   caselabel='cases n={}'.format(len(classifier.cases)),
                                   controllabel='controls n={}'.format(len(classifier.controls)),
                                   ytitle='fraction of reads', xtitle='methylation density', title='normalized read counts',
                                   legend=True)
            countsHist.plot().savefig(
                filename + '.READ-DISTR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['EfDistributionPlot'] == True:
            print(' Returning distribution of read fractions for each MD histogram: ' + path_leaf(filename) + '.EF-DISTR.png')
            combinedEFs = pd.concat(classifier.readEFsPerMDtables, axis=1)
            efsHist = histogram(combinedEFs, cases=classifier.cases_filt,
                                controls=classifier.controls_filt, meCpGs=len(classifier.CpGcolumns),
                                caselabel='cases n={}'.format(len(classifier.cases)),
                                controllabel='controls n={}'.format(len(classifier.controls)),
                                ytitle='fraction of reads', xtitle='methylation density', title='normalized read fractions',
                                legend=True)
            efsHist.plot().savefig(
                filename + '.EF-DISTR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['readcountheatmaps'] == True:
            print(' Returning TPR/FPR/TPR-FPR heatmaps using read counts: ' + path_leaf(filename) + '.COUNT-*.png')
            countTPR, countFPR, countDIFF = classifier.buildMatrices(
                read_metric='count')
            countTPRheatmap = heatmap(
                matrix=countTPR, xticks=classifier.readCountCutoffRange, yticks=classifier.mdcutoffs,
                xtitle='sample read count cutoff', ytitle='MD cutoff', title='TPR')
            countFPRheatmap = heatmap(
                countFPR, xticks=classifier.readCountCutoffRange, yticks=classifier.mdcutoffs,
                xtitle='sample read count cutoff', ytitle='MD cutoff', title='FPR')
            countDIFFheatmap = heatmap(
                countDIFF, xticks=classifier.readCountCutoffRange, yticks=classifier.mdcutoffs,
                xtitle='sample read count cutoff', ytitle='MD cutoff', title='TPR - FPR')

            countTPRheatmap.plot().savefig(
                filename + '.COUNT-TPR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            countFPRheatmap.plot().savefig(
                filename + '.COUNT-FPR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            countDIFFheatmap.plot().savefig(
                filename + '.COUNT-DIFF.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['Efheatmaps'] == True:
            print(' Returning TPR/FPR/TPR-FPR heatmaps using read fractions: ' + path_leaf(filename) + '.EF-*.png')
            efTPR, efFPR, efDIFF = classifier.buildMatrices(
                read_metric='fraction')
            efTPRheatmap = heatmap(
                efTPR, xticks=classifier.readEFCutoffRange, yticks=classifier.mdcutoffs,
                xtitle='sample read fraction cutoff', ytitle='MD cutoff', title='TPR')
            efFPRheatmap = heatmap(
                efFPR, xticks=classifier.readEFCutoffRange, yticks=classifier.mdcutoffs,
                xtitle='sample read fraction cutoff', ytitle='MD cutoff', title='FPR')
            efDIFFheatmap = heatmap(
                efDIFF, xticks=classifier.readEFCutoffRange, yticks=classifier.mdcutoffs,
                xtitle='sample read fraction cutoff', ytitle='MD cutoff', title='TPR - FPR')

            efTPRheatmap.plot().savefig(
                filename + '.EF-TPR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            efFPRheatmap.plot().savefig(
                filename + '.EF-TPR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()
            efDIFFheatmap.plot().savefig(
                filename + '.EF-TPR.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
            plt.close()

        if args_dict['optimalMDreadcounts'] == True:
            print(' Returning sample read counts above optimal MD: ' + path_leaf(filename) + '.READ-COUNTS-OPTMD.csv')
            countPerformance = classifier.cutoffPerformance(read_metric='count')
            countOptMD = classifier.optimalMDcutoff(countPerformance)
            countOptMDVals = classifier.sampleValsForMD(countOptMD)[0]
            countOptMDVals.to_csv(filename + '.READ-COUNTS-OPTMD.csv')

        if args_dict['optimalMDEf'] == True:
            print(' Returning sample read fractions above optimal MD: ' + path_leaf(filename) + '.EFS-OPTMD.csv')
            efPerformance = classifier.cutoffPerformance(read_metric='fraction')
            efOptMD = classifier.optimalMDcutoff(efPerformance)
            efOptMDVals = classifier.sampleValsForMD(efOptMD)[1]
            efOptMDVals.to_csv(filename + '.EFS-OPTMD.csv')

        print (' ')
        if outdir != None:
            print('Files stored in: {}'.format(outdir))
        else:
            print('Files stored in: {}'.format(cwd))
        print(' ')
        print('MDBC analysis completed.')
        print(' ')

if __name__ == '__main__':

	main()




