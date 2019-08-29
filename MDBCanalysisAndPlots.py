#!/usr/bin/env	python

import pandas as pd
import numpy as np
import math
import sys
import argparse
from argparse import RawTextHelpFormatter

import datetime
import ntpath

from sklearn.metrics import roc_curve, auc
import scipy.stats as stats
from scipy.stats import mannwhitneyu as ranksum

from matplotlib import pyplot as plt
import matplotlib as mpl


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Command line variables:


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
	description=
	'''
	Performs MDBC analysis on methylationDensityTable.csv.
	Computes optimal methylation density (MD) and epiallelic fraction (EF) cutoffs.
	Returns summary table of AUC, optimal EF cutoff, sensitivity, and specificity for each MD cutoff.
	Returns ROC curve and boxplots of overall optimal MD cutoff for cases vs controls.

	If indicated, can return additional plots or average methylation analysis. See flags below.
	
	Some calculations defined:
	Average methylation:
	number of meCpGs in sample reads / total CpGs in all sample reads

	EF = epiallelic fraction. Defined as:
	number of CpGs covered by sample reads with MD >= cutoff / total CpGs in all sample reads
	(special case: if MD cutoff = 0, then choose reads with MD > 0)

	if '-C' used, then EF:
	number of CpGs covered by sample reads with MD at or above MD cutoff / number of CpGs in one read.
	(NOTE: Useful mainly for DREAMing or amplicon sequencing data, in which all reads expected to cover same number of CpGs)

	The interpretation of EF in the outputs depends on the flags -C, -F sampleFractions.csv.
	For instance:
	EF can represent the sample fraction of reads at or above a given MD.
	With -C: The EF can also represent the number of reads (per sample) at or above a given MD.
	With -F: The EF can also represent the fraction of reads at or above a given MD per sampleInput.
	With -C and -F: The EF can also represent the number of reads at or above a given MD per sample input.
	For example, number of reads per mg stool.
	Or number of reads per 500mg of stool (assuming sample input the same for all samples)

	Input2Bg = False ignores the background input copies used for loading estimating number of Genomic DNA Equivalents loaded into a DREAMing assay.
	As such, EFs could be fractions of reads relative to the recorded reads in a sample and not adjusted to the estimated Genomic DNA Equivalents loaded.
	(This is an option when making the methylationDensityTable.csv)
	''',

	epilog=
	'''
	example usage:
	./MDBCanalysisAndPlots.py -table methylation_density_table.ZNF154reads.sam.csv -S MDBC_samples.csv -T Lung -N Control -P
	./MDBCanalysisAndPlots.py -table methylation_density_table.tom_dreaming.csv -T set3  -N set1 -mP -Ti .3 -Ni .1

	NOTE:
	In addition to labels for case and control sets (-T/-N), arguments require either a sampleTable.csv, and/or caseID and ctrID arguments.
	Otherwise case and control samples cannot be selected from the methylationDensityTable.csv.
	'''
	)

parser.add_argument('-table', '--methylation', required=True, metavar='methylationDensityTable.csv', action="store",
	help='''REQUIRED.
	Methylation Density Table with counts of reads at different methylation densities (MD) and corresponding number of unmethylated (numU) and methylated (numM) CpGs for each sample.
	Generated from convert[XXX]ToMethDensTable.py, where [XXX] could be 'SamReads', or 'RawDREAMing', depending on the raw data type.
	
	example of a methylationDensityTable.csv:
	numU  numM    	MD  sample1   sample2
	1     0  	0.00  9579.0  39289.0
	0     1  	1.00     6.0      0.0
	12     2  	0.14    10.0     18.0
	11     3  	0.21     2.0      5.0
	10     4  	0.28     2.0      3.0
	'''
	)

parser.add_argument('-S', '--samples', metavar='sampleTable.csv', action="store",
	help='''
	Sample Table where each column is a sample set and the corresponding rows are names of samples within the set.
	ex:
	LungPlasmaWGBS	LiverTissueWGBS
	L1		H1	
	L2		H2
	'''
	)

parser.add_argument('-T', '--cases', required=True, metavar='caseSet', action="store",
	help='''REQUIRED.
	Name of sample set to use as 'cases'.
	'''
	)

parser.add_argument('-Ti', '--caseID', metavar='caseID',
	help='''
	Identifier of case samples to use that is within the sample names in the Methylation Density Table.
	For example: If case sample names in methylationDensityTable.csv end with .3, then: -ti .3 will select them as cases.
	'''
	)

parser.add_argument('-N', '--controls', required=True, metavar='controlSet', action="store",
	help='''REQUIRED.
	Name of sample set to use as 'control'.
	'''
	)

parser.add_argument('-Ni', '--ctrID', metavar='controlID', action="store",
	help='''
	Identifier of control samples to use that is within the sample names in the Methylation Density Table.
	For example: If control sample names in methylationDensityTable.csv end with .1, then: -ti .1 will select them as cases.
	'''
	)

parser.add_argument('-F', '--fractions', metavar='sampleFractions.csv', action="store",
	help='''
	Table with relative sampleInput fractions.
	If supplied, then will perform sample adjustments based on fractions in table.
	First column is sample names ('samples') of a given sample set and second column is sample fractions ('fractions') within the set.
	ex:
	samples    fractions
	A 		1.0
	B 		0.75
	C 		0.95
	...

	OR

	samples    fractions
	A 		55.0
	B 		102.6
	C 		20.0
	...
	where the fractions now represent actual mg amounts of stool processed per sample.


	IMPORTANT: The fractions could all be based on the relative fraction of starting material used.
	As in, if 500mg of stool was used for each sample, then the fractions could be the fraction of this
	500mg used per sample. Thus sample EFs (or reads) would be normalized based on 500mg of stool.

	Alternatively, fractions could be changed to the actual amount of material used.
	For example, 50mg stool, 125mg stool, 75.5mg of stool; or 1mL plasma, 0.75mL plasma, etc.
	Then, the EFs (or reads) would be reads/sampleInput, where sampleInput could be mg stool, mLs plasma, etc.

	Could also be EF/sampleInput, too (if useCounts=False, and thus EF = fraction of reads at given MD per mg stool, or sampleInput)
	'''
	)

parser.add_argument('-C', '--useCounts', action='store_true', default=False,
	help='''
	Sets the script to use counts of reads instead of epiallelic fractions.
	'''
	)

#--------------------------------------------------------------------------
# plotting options:

#--------------------------------------------------------------------------
# visualizing the read composition of the samples:
parser.add_argument('-tR', '--totalReads', action='store_true', default=False,
	help='''
	Returns boxplot of total read EFs in cases and control samples. (All recorded reads irrespective of MD)
	This could be total reads, or total reads normalized to sampleInput.
	For example, a sample could have 0.8 total reads/mg stool.
	Or 450 reads per sampleInput (if all sample started with same amount).
	Otherwise just total reads recorded for each sample.
	-C does not affect this.

	Reports classification performance of using only the total read EF.
	'''
	)

parser.add_argument('-rF', '--readEFperMD', action='store_true', default=False,
	help='''
	Returns image showing cases vs controls read EFs for each MD. (EFs = fraction of reads with exactly a specific MD)
	This is the fraction of reads (EF) with a given MD out of the total reads for each sample.
	Could be normalized to sampleInput.
	For example, a sample could have an EF of reads with MD=0.2 of 0.8 per mg stool.
	Or 450 reads with MD=0.2 per sampleInput (if all sample started with same amount).
	-C does not affect this.

	Uses MDs of reads recorded in the Methylation Density Table.
	'''
	)

parser.add_argument('-rP', '--readsPerMDPlot', action='store_true', default=False,
	help='''
	Returns stacked bar plot showing read counts per input amount.
	read counts could be sample counts of reads at each MD, or sample counts of reads at each MD normalized to sampleInput.
	SampleInput amount based on the amount of sample analyzed.
	Most relevant if -F 'sampleFractions.csv' given, where the fractions are actual amounts of sample assessed (mLs plasma, mgs of stool, etc).
	For example, a sample could have 0.8 reads with MD = 1.0 per mg stool.
	Or 450 reads with MD = 1.0 per sampleInput (if all sample started with same amount)
	Otherwise just recorded read counts at each MD is plotted per sample.
	-C does not affect this.
	'''
	)

parser.add_argument('-wP', '--weightedReadFractionPlot', action='store_true', default=False,
	help='''
	Returns stacked bar plot of weighted locus methylation density read fractions for each sample.

	Weighted locus methylation density fraction =
	number of C's covered by reads with given methylation density / total C's covered by reads in sample

	NOTE that the profile is based on the relative fraction of reads within a sample. Not compared between samples.
	So adjustment by the -F sampleFractions.csv will not change the profile, as the relative fraction of reads in each sample will still stay the same.
	However, the Input2Bg (background reads included as reads with MD=0; option when making the methylationDensityTable.csv) will change the relative fraction of reads in a sample and thus the profile of this plot.
	-C does not affect this.

	Samples ordered by average locus methylation (total number of meC's covered by reads / total C's covered by reads).
	'''
	)

parser.add_argument('-his', '--readEFhistogram', action='store_true', default=False,
	help='''
	Returns histogram of pooled sample read EFs across MDs for each sample set. Does not show total reads, only relative read EFs within a sample set.

	NOTE that the profile is based on the relative fraction of reads for all pooled samples in a given sample set.
	Adjustment by the -F sampleFractions.csv may change the profile if the relative fractions of samples vary considerably between each other.
	-C does not affect this.
	'''
	)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# the MDBC analysis:

parser.add_argument('-mP', '--MDBCplots', action='store_true', default=False,
	help='''
	Returns ROC curves and boxplots for each MD cutoff during the methylation density cutoff iteration, not just the optimal.

	-C and -F do affect this.
	Because using counts of reads at or above a given MD, or the sample fraction of these reads, and adjusting these counts or fractions to the amount of sampleInput can all change the classification performance.

	Again, the EF can be interpreted as the count of reads, or fraction of reads, at or above a given MD, normalized to the sampleInput fraction depending on the -C and -F flags. 
	'''
	)

parser.add_argument('-heat', '--heatmaps', action='store_true', default=False,
	help='''
	Returns TPR, FPR, and TPR-FPR (Diff) heatmaps for each iterative MD/EF cutoff combination.
	'''
	)

parser.add_argument('-aP', '--averageMethylation', action='store_true', default=False,
	help='''
	Include average methylation analysis. Plots average methylation ROC and returns average methylation boxplot (cases vs. controls).
	'''
	)
#--------------------------------------------------------------------------

parser.add_argument('-out', '--outDir', action='store',
	help='''
	Directory path to store output files.
	'''
	)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

#aesthetics for plots:
font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 12}
mpl.rc('font', **font)


#--------------------------------------------------------------------------

# get file name from path no matter the operating system or file path
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Set up global variables:


# Collect command line arguments
args = parser.parse_args()

# Date run:
date = datetime.date.today()


# REQUIRED ARGS:
# methylationDensityTable.csv
df = pd.read_csv(args.methylation)
if df.empty == True:
	print 'The loaded methylationDensityTable.csv is empty.'
	print 'check the -M input file.'
	print "Exiting script."
	sys.exit(0)

# Name of methylationDensityTable.csv
inputName = path_leaf(str(args.methylation))

# locus of interest (should be contained in methylationDensityTable.csv name)
locus = str(inputName).split('.', 1)[1].split('.Created')[0]
# methylation_density_table.tom_dreaming.Input2Bg=True.Created=2019-08-07.csv --> tom_dreaming.Input2Bg=True


# cases sample set name
if (vars(args)['cases'] == None):
	print 'Need to indicate a name for the cases sample set with -T cases.'
	print "Exiting script."
	sys.exit(0)
else:
	casesType = str(args.cases)


# controls sample set name
if (vars(args)['controls'] == None):
	print 'Need to indicate a name for the controls sample set with -T controls.'
	print "Exiting script."
	sys.exit(0)
else:
	controlsType = str(args.controls)


#--------------------------------------------------------------------------
# Select samples:


# No sampleTable.csv and no caseID or ctrID given, then cannot select samples and exit script
if (vars(args)['samples'] == None) and ((vars(args)['ctrID'] == None) or (vars(args)['caseID'] == None)):
	print '''Need to select case and control groups either using 'sampleTable.csv' or caseID/ctrID'''
	sys.exit(0)


# if sampleTable.csv is available:
if vars(args)['samples'] != None:

	sampleTable = pd.read_csv(args.samples)
	sampleSets = sampleTable.columns.values.tolist()

	# check if casesType in sampleTable.csv
	if casesType in sampleSets:
		cases = [str(i) for i in sampleTable[casesType].values if str(i) != 'nan' ]

		# if caseIDs also exist to futher select samples...
		if vars(args)['caseID'] != None:
			ID = vars(args)['caseID']
			cases = [i for i in cases if str(ID) in i]

	# if not, then presumably select cases just using caseID from columns in methylationDensityTable.csv
	elif vars(args)['caseID'] != None:
		ID = vars(args)['caseID']
		cases = df.loc[:, [col for col in df.columns if str(ID) in col] ].columns.values

	else:
		print str(casesType) + " not found in: " + str(args.samples)
		print "and caseID not given. Cannot select cases."
		print "Exiting script."
		sys.exit(0)


	# check if controlsType in sampleTable.csv
	if controlsType in sampleSets:
		controls = [str(i) for i in sampleTable[controlsType].values if str(i) != 'nan' ]

		# if ctrIDs also exist to futher select samples...
		if vars(args)['ctrID'] != None:
			ID = vars(args)['ctrID']
			controls = [i for i in controls if str(ID) in i]

	# if not, then presumably select controls just using ctrID from columns in methylationDensityTable.csv
	elif vars(args)['ctrID'] != None:
		ID = vars(args)['ctrID']
		controls = df.loc[:, [col for col in df.columns if str(ID) in col] ].columns.values

	else:
		print str(controlsType) + " not found in: " + str(args.samples)
		print "And ctrID also not given. Cannot select cases."
		print "Exiting script"
		sys.exit(0)


# if case or ctr IDs are indicated but no sampleTable.csv:
else:
	
	if vars(args)['ctrID'] != None:
		ID = vars(args)['ctrID']
		controls = df.loc[:, [col for col in df.columns if str(ID) in col] ].columns.values

	if vars(args)['caseID'] != None:
		ID = vars(args)['caseID']
		cases = df.loc[:, [col for col in df.columns if str(ID) in col] ].columns.values

#--------------------------------------------------------------------------

# --fractions ; sampleFractions.csv
inputFracAdjust = False
if vars(args)['fractions'] != None:
	inputFracs_df = pd.read_csv(vars(args)['fractions'])
	inputFracs = dict(zip( [str(i) for i in inputFracs_df['samples'].values], inputFracs_df['fractions'].values))
	inputFracAdjust = True

# '-C', '--useCounts' flag
usePeakCounts = False
if vars(args)['useCounts'] == True:
	usePeakCounts = True
	numCpGs = (df['numU'] + df['numM']).mean()


if vars(args)['outDir'] != None:
	outdir = str(args.outDir)
	fileTag = outdir + 'MDBCAnalysis.' + 'Created=' + str(date) + '.' + casesType + '_vs_' + controlsType + '.' + locus  + '.inputFracAdjust=' + str(inputFracAdjust) + '.usePeakCounts=' + str(usePeakCounts)
else:
	# ex: MDBCAnalysis.Created=2019-08-09.set3_vs_set1.tom_dreaming.Input2Bg=False.inputFracAdjust=True.usePeakCounts=True.FILEtype
	fileTag = 'MDBCAnalysis.' + 'Created=' + str(date) + '.' + casesType + '_vs_' + controlsType + '.' + locus  + '.inputFracAdjust=' + str(inputFracAdjust) + '.usePeakCounts=' + str(usePeakCounts)


# print std.out 'print' statements to new file and to console:
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(fileTag + '.out', "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        pass

sys.stdout = Logger()


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

print 'Date:' + str(date)
print '#--------------------------------------------------------------------------'
print 'Command:'
print sys.argv
print ' '
print args
print ' '
print '#--------------------------------------------------------------------------'
print 'MDBC ANALYSIS BEGIN'
print 'Methylation Density Table: ' + str(inputName)
print 'inputFracAdjust = ' + str(inputFracAdjust)
if inputFracAdjust == True:
	print 'inputFracAdjust file: ' + vars(args)['fractions']
print 'usePeakCounts = ' + str(usePeakCounts)
print 'Cases = ' + casesType
print 'Controls = ' + controlsType
print ' '


# lists to store AUC, optimal EF, and associated sensitivity and specificity for each MD
aucs_per_md = []
opt_df_thresh_per_md = []
sens_per_md = []
spec_per_md = []


# lists to build TPR and FPR matrices for plotting heat maps
tpr_matrix = []
fpr_matrix = []


# methylation density (md) cutoffs (0%, 5%, 10%, etc...) to 100%:
md_cutoffs = np.arange(0,1.05,0.05)


#--------------------------------------------------------------------------
# Set EF cutoff range for heat maps and MDBC based on EFs (or potentially counts of methylated reads) in sample sets
# Also remove samples with no reads


# Counts of methylated reads and total reads in samples:
methReadCounts = df[ ['MD'] + list(cases) + list(controls) ].copy().set_index('MD').loc[ [i for i in set(df['MD'].values) if i != 0 ] ].sum()
TotalReadCounts = df[ list(cases) + list(controls) ].copy().sum()


# remove any samples that have no reads at all:
noValueSamples = list(TotalReadCounts[TotalReadCounts == 0].index)

print 'Samples with no reads: '
print noValueSamples
print 'Removed'
print ' '

methReadCounts_filt = methReadCounts[ ~methReadCounts.index.isin(noValueSamples) ].copy()
TotalReadCounts_filt = TotalReadCounts[ ~TotalReadCounts.index.isin(noValueSamples) ].copy()

cases = [i for i in cases if i not in noValueSamples ]
controls = [i for i in controls if i not in noValueSamples ]

casesAndControls = list(cases) + list(controls)

if inputFracAdjust == True:
	sampleFracs = [float(inputFracs[i]) for i in casesAndControls]
	casesFracs = [float(inputFracs[i]) for i in cases]
	ctrlFracs = [float(inputFracs[i]) for i in controls]
else:
	sampleFracs = [1.0 for i in casesAndControls]
	casesFracs = [1.0 for i in cases]
	ctrlFracs = [1.0 for i in controls]

#--------------------------------------------------------------------------
# set epifraction (ef) cutoffs
if usePeakCounts == True:
	methReadEFs = methReadCounts_filt 
	# note that for MD cutoff of 0, reads of MD > 0 used, so only methylated reads are of interest.
	# Therefore, just using total reads (which would include reads with MD = 0) not needed to set thresh
	# and also samples not classified using just total amount of reads recorded

else:
	methReadEFs = methReadCounts_filt / TotalReadCounts_filt


#--------------------------------------------------------------------------
# Set EF cutoff limit for plotting on heatmaps

if inputFracAdjust == True:

	# if true, adjust the EFs
	methReadEFs = methReadEFs / sampleFracs

	# EFs could be fracs or positive numbers depending on the fraction adjustment
	max_ef = max(methReadEFs)*1.10

else:

	# if using peak counts, then EF cutoff will be positive integers above one, so adjust accordingly
	if usePeakCounts == True:
		max_ef = max(methReadEFs)*1.10
	else:
		# otherwise, EFs will be fractions, and limit to 1.0
		if max(methReadEFs)*1.10 > 1.0:
			max_ef = 1.0
		else:
			max_ef = max(methReadEFs)*1.10


ef_cutoffs = np.array(list(np.linspace(0.000, max_ef, 100, endpoint=False)) + [max_ef])


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' Sample summary: overall sample read fractions per MD and total reads per sample set '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Total Reads in Sample Sets

if (vars(args)['totalReads'] == True):

	print 'Generating total read EFs per sample set...'


	if inputFracAdjust == True:
		print 'Adjusted by relative sample input fractions'
		a = list(TotalReadCounts_filt[TotalReadCounts_filt.index.isin(cases)].copy().values / casesFracs)
		b = list(TotalReadCounts_filt[TotalReadCounts_filt.index.isin(controls)].copy().values / ctrlFracs)
	else:
		a = TotalReadCounts_filt[TotalReadCounts_filt.index.isin(cases)].copy().values.tolist()
		b = TotalReadCounts_filt[TotalReadCounts_filt.index.isin(controls)].copy().values.tolist()


	fig, ax = plt.subplots(figsize=(2,4))
	bp = ax.boxplot([a, b], positions=[1,3], widths=1, patch_artist=True)

	for flier in bp['fliers']:
		flier.set(marker='', color='black')
	for whisker in bp['whiskers']:
		whisker.set(color='black', linewidth=1)
	for cap in bp['caps']:
		cap.set(color='black', linewidth=1)
	for median in bp['medians']:
		median.set(color='black', linewidth=1)

	bp['boxes'][0].set( color='red', facecolor='salmon', linewidth=2, alpha=0.5)
	bp['boxes'][1].set( color='blue', facecolor='skyblue', linewidth=2, alpha=0.5)

	scatter = ax.scatter(x = np.random.normal(1, 0.1, size=len(a)),
		y = a, c='r', marker='.', edgecolors='', s = 50)
	scatter = ax.scatter(x = np.random.normal(3, 0.1, size=len(b)),
		y = b, c='b', marker='.', edgecolors='', s = 50)

	s, p = stats.ranksums(a, b)
	print 'Boxplots Case vs Control Total Read EFs:'
	print '	# of cases = ' + str(len(a))
	print '	# of controls = ' + str(len(b))
	print '	Case median total read EF = ' + str(np.median(a))
	print '	Control median total read EF = ' + str(np.median(b))
	print '	Case lowest total read EF = ' + str(np.min(a))
	print '	Control lowest total read EF = ' + str(np.min(b))
	print '	Case highest total read EF = ' + str(np.max(a))
	print '	Control highest total dread EF = ' + str(np.max(b))
	print '	ranksum p-val = ' + str(p)

	# Identify optimal thresholds, plot, and print sens/spec:
	case_labels = [1 for i in a]
	ctrl_labels = [0 for i in b]
	roc_values = a + b
	roc_labels = case_labels + ctrl_labels
	fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
	roc_auc = auc(fpr, tpr)
	optimal_idx = np.argmax( tpr - fpr )
	optimal_threshold = thresholds[optimal_idx]
	#print '	optimal total read threshold = ' + str(optimal_threshold)
	print '	Classification by total read EF:'
	print '	AUC = ' + str(round(roc_auc, 2))
	print '	TPR = ' + str(round(tpr[optimal_idx], 2))
	print '	1-FPR = ' + str(round(1.0-fpr[optimal_idx], 2))
	print '	optimal total read threshold = ' + str(optimal_threshold)
	print ' '


	#ax.axhline(y=optimal_threshold, linestyle='--', color='k')

	max_val = np.max(roc_values)

	s, p = stats.ranksums(a, b)
	if 0.01 <= p < 0.05:
		plt.text(x=1.7,y=max_val * 1.10, s='*', fontsize=30)
	if 0.001 <= p < 0.01:
		plt.text(x=1.6,y=max_val * 1.10, s='**', fontsize=30)
	if p < 0.001:
		plt.text(x=1.5,y=max_val * 1.10, s='***', fontsize=30)
	if p >= 0.05:
		plt.text(x=1.5,y=max_val * 1.12, s='ns', fontsize=30)


	ax.set_xlim([0, 4])
	ax.set_ylabel('total read EF', fontsize=18, rotation=0, labelpad=40)
	ax.yaxis.set_label_coords(-0.97,0.43)
	ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
	ax.set_ylim([0, max_val * 1.25])
	plt.title('sample reads', fontsize=20)
	plt.yticks(np.linspace(0, max_val * 1.15, 5), ['0.0'] + [ '%.1e' % i for i in np.linspace(0, max_val * 1.15, 5)[1:] ], fontsize=18)
	fig.savefig(fileTag + '.TOTALREAD' + '_BOXPLOT.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
	plt.close()

#--------------------------------------------------------------------------

# Sample read fractions per MD cutoff Boxplot

if (vars(args)['readEFperMD'] == True):

	print 'Calculating sample read EFs for each MD ...'
	if inputFracAdjust == True:
		print 'Adjusted by relative sample input fractions'

	sorted_MDs = sorted(list(set(df['MD'].values)))

	ncols = 10
	nrows = math.ceil(len(sorted_MDs) / float(ncols))

	fig = plt.figure( figsize=(5*ncols, 4*nrows) )

	x = 1
	for md in sorted_MDs:

		md_counts = df[ ['MD'] + casesAndControls ].copy().set_index('MD').loc[md]

		# several rows could have been selected if MD index happened to be the same
		# If so, sum counts in rows together into a series object. Otherwise, keep the series
		if len(md_counts.shape) > 1:
			sample_read_counts_at_md = md_counts.sum()
		else:
			sample_read_counts_at_md = md_counts

		# if no reads at given MD and thus NaN values, replace with 0
		sample_read_counts_at_md.fillna(0, inplace=True)

		if inputFracAdjust == True:
			EF_at_MD = (sample_read_counts_at_md / TotalReadCounts_filt) / sampleFracs
		else:
			EF_at_MD = (sample_read_counts_at_md / TotalReadCounts_filt)


		a = EF_at_MD[EF_at_MD.index.isin(cases)].copy().values.tolist()
		b = EF_at_MD[EF_at_MD.index.isin(controls)].copy().values.tolist()

		ax = fig.add_subplot(nrows, ncols, x)

		bp = ax.boxplot([a, b], positions=[1,3], widths=1, patch_artist=True)

		for flier in bp['fliers']:
			flier.set(marker='', color='black')
		for whisker in bp['whiskers']:
			whisker.set(color='black', linewidth=1)
		for cap in bp['caps']:
			cap.set(color='black', linewidth=1)
		for median in bp['medians']:
			median.set(color='black', linewidth=1)

		bp['boxes'][0].set( color='red', facecolor='salmon', linewidth=2, alpha=0.5)
		bp['boxes'][1].set( color='blue', facecolor='skyblue', linewidth=2, alpha=0.5)

		scatter = ax.scatter(x = np.random.normal(1, 0.1, size=len(a)),
			y = a, c='r', marker='.', edgecolors='', s = 50)
		scatter = ax.scatter(x = np.random.normal(3, 0.1, size=len(b)),
			y = b, c='b', marker='.', edgecolors='', s = 50)

		s, p = stats.ranksums(a, b)
		print ' '
		print '	Only reads with MD = ' + str(md)
		print '		# of cases = ' + str(len(a))
		print '		# of controls = ' + str(len(b))
		print '		Case median EF = ' + str(np.median(a))
		print '		Control median EF = ' + str(np.median(b))
		print '		Case lowest EF = ' + str(np.min(a))
		print '		Control lowest EF = ' + str(np.min(b))
		print '		Case highest EF = ' + str(np.max(a))
		print '		Control highest EF = ' + str(np.max(b))
		print '		ranksum p-val = ' + str(p)

		# Identify optimal thresholds, plot, and print sens/spec:
		case_labels = [1 for i in a]
		ctrl_labels = [0 for i in b]
		roc_values = a + b
		roc_labels = case_labels + ctrl_labels
		fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
		roc_auc = auc(fpr, tpr)
		optimal_idx = np.argmax( tpr - fpr )
		optimal_threshold = thresholds[optimal_idx]
		print '		optimal EF threshold = ' + str(optimal_threshold)

		#ax.axhline(y=optimal_threshold, linestyle='--', color='k')

		max_val = float(np.max(roc_values))

		if max_val != 0.0:

			s, p = stats.ranksums(a, b)
			if 0.01 <= p < 0.05:
				ax.text(x=1.7,y=max_val * 1.10, s='*', fontsize=30)
			if 0.001 <= p < 0.01:
				ax.text(x=1.6,y=max_val * 1.10, s='**', fontsize=30)
			if p < 0.001:
				ax.text(x=1.5,y=max_val * 1.10, s='***', fontsize=30)
			if p >= 0.05:
				ax.text(x=1.5,y=max_val * 1.12, s='ns', fontsize=30)

			ax.set_xlim([0, 4])
			ax.set_ylabel('epiallelic fraction', fontsize=20)
			#ax.yaxis.set_label_coords(-0.80,0.43)
			ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
			ax.set_ylim([0, max_val * 1.25])
			ax.set_title('MD = ' + str(round(md, 3)), fontsize=20)
			ax.set_yticks(np.linspace(0, max_val * 1.15, 5))
			ax.set_yticklabels(['0.0'] + [ '%.1e' % i for i in np.linspace(0, max_val * 1.15, 5)[1:] ], fontsize=18)

		else: # If no reads at MD and thus EF = 0
			max_val = 0.0001

			s, p = stats.ranksums(a, b)
			if 0.01 <= p < 0.05:
				ax.text(x=1.7,y=max_val * 1.10, s='*', fontsize=30)
			if 0.001 <= p < 0.01:
				ax.text(x=1.6,y=max_val * 1.10, s='**', fontsize=30)
			if p < 0.001:
				ax.text(x=1.5,y=max_val * 1.10, s='***', fontsize=30)
			if p >= 0.05:
				ax.text(x=1.5,y=max_val * 1.12, s='ns', fontsize=30)

			ax.set_xlim([0, 4])
			ax.set_ylabel('epiallelic fraction', fontsize=20)
			#ax.yaxis.set_label_coords(-0.80,0.43)
			ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
			ax.set_ylim([0, max_val * 1.25])
			ax.set_title('MD = ' + str(round(md, 3)), fontsize=20)
			ax.set_yticks(np.linspace(0, max_val * 1.15, 5))
			ax.set_yticklabels(['0.0'] + [ '%.1e' % i for i in np.linspace(0, max_val * 1.15, 5)[1:] ], fontsize=18)

		x += 1

	print ' '
	print 'Generating sample set read fraction per MD Boxplot image...'
	print str(len(sorted_MDs)) + ' MD boxplots to plot...'

	plt.tight_layout()
	fig.savefig(fileTag + '.ReadFracPerMD-BOXPLOT.png', bbox_inches='tight', pad_inches=1.0, dpi=400)
	plt.close()


	print 'Sample read EF per MD boxplot completed'
	print ' '


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' Stacked Barplots of read MD per input '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

if (vars(args)['readsPerMDPlot'] == True):

	print 'Generating stacked bar plots of read counts per MD for each sample set...'

	if inputFracAdjust == True:
		print 'Adjusted by relative sample input fraction'


	# list of observed read methylation densities that exist in the methylation density table
	sorted_MDs = sorted(list(set(df['MD'].values)))
	sorted_MDs = [ i for i in sorted_MDs if i != 0.0 ] # don't care about unmethylated reads, only care about signal coming from methylated

	# table with meC and total C information for each read in table
	MD_cols = df[['numU', 'numM', 'MD']].copy()
	MD_cols['C'] = MD_cols[['numU', 'numM']].sum(axis=1)


	print 'Making sample read MD per input table for cases...'
	# construct table of read counts and their meC information for samples of interest
	cases_cols = df[cases]
	if inputFracAdjust == True:
		# adjust peak counts by sample input fraction
		cases_cols = cases_cols.copy().div(casesFracs)
	case_df = pd.concat([MD_cols, cases_cols], axis=1)
	# make new dataframe to populate with weighted sample fractions of each methylation density
	case_weighted_MDs = pd.DataFrame(index=sorted_MDs)

	# collect sample average methylation values
	caseAveMeth = []

	for case in cases:
		temp = MD_cols.copy() # temporary table of meC and total C information of read patterns
		counts = case_df[case] # counts of reads for each different pattern in case

		# get total number of C's covered by reads in case:
		temp['sample_C'] = temp['C'].values * counts.values
		total_C = np.nansum(temp['sample_C'])

		# compute average methylation (use to sort samples in plot)
		# (total number of meC's covered by reads / total C's covered by reads)
		ave_meC = np.nansum(temp['numM'].values * counts.values) / float(total_C)
		caseAveMeth.append(ave_meC)

		counts = [] # populate list with case fraction of reads with a given methylation density
		for m in sorted_MDs:
			selection = case_df[case_df['MD'] == m][case].values # get counts of reads with given MD for the sample

			counts.append(sum(selection))

		case_weighted_MDs[case] = counts # append fractions to dataframe for given case

	#--------------------------------------------------------------------------
	# reorder samples (columns) based on average methylation (meC of sample)
	reorder = zip(caseAveMeth, cases)
	reorder = sorted(reorder, key = lambda x: x[0])
	reorder_names = [i[1] for i in reorder]
	case_weighted_MDs = case_weighted_MDs[reorder_names]
	#--------------------------------------------------------------------------

	print 'Making weighted locus methylation density table for controls...'
	ctrl_cols = df[controls]
	if inputFracAdjust == True:
		ctrl_cols = ctrl_cols.copy().div(ctrlFracs)
	ctrl_df = pd.concat([MD_cols, ctrl_cols], axis=1)
	# make new dataframe to populate with weighted sample fractions of each methylation density
	ctrl_weighted_MDs = pd.DataFrame(index=sorted_MDs)

	# collect sample average methylation values
	ctrlAveMeth = []

	for ctrl in controls:
		temp = MD_cols.copy() # temporary table of meC and total C information of read patterns
		counts = ctrl_df[ctrl] # counts of reads for each different pattern in ctrl

		# get total number of C's covered by reads in ctrl:
		temp['sample_C'] = temp['C'].values * counts.values
		total_C = np.nansum(temp['sample_C'])

		# compute average methylation
		# (total number of meC's covered by reads / total C's covered by reads)
		ave_meC = np.nansum(temp['numM'].values * counts.values) / float(total_C)
		ctrlAveMeth.append(ave_meC)

		counts = [] # populate list with control fraction of reads with a given methylation density
		for m in sorted_MDs:
			selection = ctrl_df[ctrl_df['MD'] == m][ctrl].values # get counts of reads with given MD for the sample

			counts.append(sum(selection))

		ctrl_weighted_MDs[ctrl] = counts # append fractions to dataframe for given ctrl

	#--------------------------------------------------------------------------
	# reorder samples (columns) based on average methylation (meC of sample)
	reorder = zip(ctrlAveMeth, controls)
	reorder = sorted(reorder, key = lambda x: x[0])
	reorder_names = [i[1] for i in reorder]
	ctrl_weighted_MDs = ctrl_weighted_MDs[reorder_names]
	#--------------------------------------------------------------------------
	#--------------------------------------------------------------------------

	print 'Plotting weighted locus methylation density for cases...'
	# Plotting
	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])

	cmap = plt.get_cmap('coolwarm')
	colors = cmap(sorted_MDs) # colors based on the MDs, which should be between 0 and 1.0
	plot = case_weighted_MDs.T.plot(kind='bar', stacked=True, ax=ax, color=colors, width=1, legend=True)

	plot.set_xlim( [-0.5, len(cases)-0.5] )
	ax.set_ylim(0, max(list(case_weighted_MDs.sum()) + list(case_weighted_MDs.sum()))*1.05 )

	#ax.set_yticks(np.arange(0,1.1,0.2))
	plt.yticks(fontsize=16)
	ax.set_ylabel('reads per sample input', fontsize=26)
	#ax.yaxis.set_label_coords(-0.10, 0.5)

	# make bottom x-axis labels the number of reads
	#ax.set_xticklabels([str(int(i)) for i in case_df[reorder_names].sum().values], rotation=45, fontsize=10) 
	#ax.set_xlabel('reads covering locus', fontsize=14)

	# colorbar methylation density legend
	# keep color scale from 0 to 1.0
	ax2 = fig.add_axes([0.1, 1.25, 0.5, 0.05])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1))
	cb.set_label('methylation density')
	cb.ax.xaxis.set_ticks_position('top')
	cb.set_clim(0.0, 1)

	# top x-axis labels = average methylation values
	if (vars(args)['averageMethylation'] == True):
		ax3 = ax.twiny()
		ax3.set_xticks(ax.get_xticks())
		ax3.set_xticklabels([round(i[0],3) for i in reorder], fontsize=10, rotation=45) 
		ax3.xaxis.set_ticks_position('top')
		ax3.xaxis.set_label_position('top')
		ax3.set_xlabel('average methylation', fontsize=14)
		ax3.xaxis.set_label_coords(0.5, 1.1)
		ax3.set_xlim(ax.get_xlim())

	# bottom label = sample type and number of samples
	ax4 = ax.twiny()
	ax4.set_xticks([])
	ax4.set_xticklabels([])
	ax4.set_xlabel(str(casesType) + ' (n=' + str(len(cases)) + ')', fontsize=26)
	ax4.xaxis.set_label_coords(0.5, -0.25)

	ax.get_legend().remove()

	fig.savefig(fileTag + '.CASES_FRAC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
	print 'Cases Barplot created.'
	plt.close()

	#--------------------------------------------------------------------------

	print 'Plotting weighted locus methylation density for controls...'
	# Plotting
	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])

	cmap = plt.get_cmap('coolwarm')
	colors = cmap(sorted_MDs) # colors based on the MDs, which should be between 0 and 1.0
	plot = ctrl_weighted_MDs.T.plot(kind='bar', stacked=True, ax=ax, color=colors, width=1, legend=True)

	plot.set_xlim( [-0.5, len(cases)-0.5] )
	ax.set_ylim(0, max(list(case_weighted_MDs.sum()) + list(case_weighted_MDs.sum()))*1.05 )

	#ax.set_yticks(np.arange(0,1.1,0.2))
	plt.yticks(fontsize=16)
	ax.set_ylabel('reads per sample input', fontsize=26)
	#ax.yaxis.set_label_coords(-0.10, 0.5)

	# make bottom x-axis labels the number of reads
	#ax.set_xticklabels([str(int(i)) for i in ctrl_df[reorder_names].sum().values], rotation=45, fontsize=10) 
	#ax.set_xlabel('reads covering locus', fontsize=14)

	# colorbar methylation density legend
	# keep color scale from 0 to 1.0
	ax2 = fig.add_axes([0.1, 1.25, 0.5, 0.05])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1))
	cb.set_label('methylation density')
	cb.ax.xaxis.set_ticks_position('top')
	cb.set_clim(0.0, 1)

	# top x-axis labels = average methylation values
	if (vars(args)['averageMethylation'] == True):
		ax3 = ax.twiny()
		ax3.set_xticks(ax.get_xticks())
		ax3.set_xticklabels([round(i[0],3) for i in reorder], fontsize=10, rotation=45) 
		ax3.xaxis.set_ticks_position('top')
		ax3.xaxis.set_label_position('top')
		ax3.set_xlabel('average methylation', fontsize=14)
		ax3.xaxis.set_label_coords(0.5, 1.1)
		ax3.set_xlim(ax.get_xlim())

	# bottom label = sample type and number of samples
	ax4 = ax.twiny()
	ax4.set_xticks([])
	ax4.set_xticklabels([])
	ax4.set_xlabel(str(controlsType) + ' (n=' + str(len(controls)) + ')', fontsize=26)
	ax4.xaxis.set_label_coords(0.5, -0.25)

	ax.get_legend().remove()

	fig.savefig(fileTag + '.CTRLS_FRAC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)

	print 'Control Barplot created.'
	print 'Done'
	print ' '
	plt.close()


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' Stacked Barplot of weighted sample read EFs '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

if (vars(args)['weightedReadFractionPlot'] == True):

	print 'Generating stacked bar plot of weighted locus methylation density fractions for each sample set...'

	if inputFracAdjust == True:
		print 'Adjusted by relative sample input fraction'


	# list of observed read methylation densities that exist in the methylation density table
	sorted_MDs = sorted(list(set(df['MD'].values)))

	# table with meC and total C information for each read in table
	MD_cols = df[['numU', 'numM', 'MD']].copy()
	MD_cols['C'] = MD_cols[['numU', 'numM']].sum(axis=1)


	print 'Making weighted locus methylation density table for cases...'
	# construct table of read counts and their meC information for samples of interest
	cases_cols = df[cases]
	if inputFracAdjust == True:
		# adjust peak counts by sample input fraction
		cases_cols = cases_cols.copy().div(casesFracs)
	case_df = pd.concat([MD_cols, cases_cols], axis=1)
	# make new dataframe to populate with weighted sample fractions of each methylation density
	case_weighted_MDs = pd.DataFrame(index=sorted_MDs)

	# collect sample average methylation values
	caseAveMeth = []

	for case in cases:
		temp = MD_cols.copy() # temporary table of meC and total C information of read patterns
		counts = case_df[case] # counts of reads for each different pattern in case

		# get total number of C's covered by reads in case:
		temp['sample_C'] = temp['C'].values * counts.values
		total_C = np.nansum(temp['sample_C'])

		# compute average methylation
		# (total number of meC's covered by reads / total C's covered by reads)
		ave_meC = np.nansum(temp['numM'].values * counts.values) / float(total_C)
		caseAveMeth.append(ave_meC)

		# weighted fraction
		# (number of C's covered by reads with given methylation density / total C's covered by reads in case)
		weighted_fracs = [] # populate list with case fraction of reads with a given methylation density
		for m in sorted_MDs:
			selection = temp[temp['MD'] == m] # get rows (read patterns) that have the given methylation density
			frac = np.nansum(selection['sample_C']) / float(total_C)
			weighted_fracs.append(frac)

		case_weighted_MDs[case] = weighted_fracs # append fractions to dataframe for given case

	#--------------------------------------------------------------------------
	# reorder samples (columns) based on average methylation (meC of sample)
	reorder = zip(caseAveMeth, cases)
	reorder = sorted(reorder, key = lambda x: x[0])
	reorder_names = [i[1] for i in reorder]
	case_weighted_MDs = case_weighted_MDs[reorder_names]
	#--------------------------------------------------------------------------
	print 'Plotting weighted locus methylation density for cases...'
	# Plotting
	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])

	cmap = plt.get_cmap('coolwarm')
	colors = cmap(sorted_MDs) # colors based on the MDs, which should be between 0 and 1.0
	plot = case_weighted_MDs.T.plot(kind='bar', stacked=True, ax=ax, color=colors, width=1, legend=True)

	plot.set_xlim( [-0.5, len(cases)-0.5] )
	ax.set_ylim(0,1)

	ax.set_yticks(np.arange(0,1.1,0.2))
	plt.yticks(fontsize=16)
	ax.set_ylabel('weighted sample fraction', fontsize=26)
	ax.yaxis.set_label_coords(-0.10, 0.5)

	# make bottom x-axis labels the number of reads
	ax.set_xticklabels([str('%.1e' % i) for i in case_df[reorder_names].sum().values], rotation=90, fontsize=10) 
	ax.set_xlabel('total reads per sample input', fontsize=12)

	# colorbar methylation density legend
	# keep color scale from 0 to 1.0
	ax2 = fig.add_axes([0.1, 1.25, 0.5, 0.05])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1))
	cb.set_label('locus methylation density')
	cb.ax.xaxis.set_ticks_position('top')
	cb.set_clim(0.0, 1)

	# top x-axis labels = average methylation values
	if (vars(args)['averageMethylation'] == True):
		ax3 = ax.twiny()
		ax3.set_xticks(ax.get_xticks())
		ax3.set_xticklabels([round(i[0],3) for i in reorder], fontsize=10, rotation=45) 
		ax3.xaxis.set_ticks_position('top')
		ax3.xaxis.set_label_position('top')
		ax3.set_xlabel('average methylation', fontsize=12)
		ax3.xaxis.set_label_coords(0.5, 1.1)
		ax3.set_xlim(ax.get_xlim())

	# bottom label = sample type and number of samples
	ax4 = ax.twiny()
	ax4.set_xticks([])
	ax4.set_xticklabels([])
	ax4.set_xlabel(str(casesType) + ' (n=' + str(len(cases)) + ')', fontsize=26)
	ax4.xaxis.set_label_coords(0.5, -0.25)

	ax.get_legend().remove()

	fig.savefig(fileTag + '.CASES_WEIGHTED_FRAC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
	print 'Cases Barplot created.'



	print 'Making weighted locus methylation density table for controls...'
	ctrl_cols = df[controls]
	if inputFracAdjust == True:
		ctrl_cols = ctrl_cols.copy().div(ctrlFracs)
	ctrl_df = pd.concat([MD_cols, ctrl_cols], axis=1)
	# make new dataframe to populate with weighted sample fractions of each methylation density
	ctrl_weighted_MDs = pd.DataFrame(index=sorted_MDs)

	# collect sample average methylation values
	ctrlAveMeth = []

	for ctrl in controls:
		temp = MD_cols.copy() # temporary table of meC and total C information of read patterns
		counts = ctrl_df[ctrl] # counts of reads for each different pattern in ctrl

		# get total number of C's covered by reads in ctrl:
		temp['sample_C'] = temp['C'].values * counts.values
		total_C = np.nansum(temp['sample_C'])

		# compute average methylation
		# (total number of meC's covered by reads / total C's covered by reads)
		ave_meC = np.nansum(temp['numM'].values * counts.values) / float(total_C)
		ctrlAveMeth.append(ave_meC)

		# weighted fraction
		# (number of C's covered by reads with given methylation density / total C's covered by reads in ctrl)
		weighted_fracs = [] # populate list with ctrl fraction of reads with a given methylation density
		for m in sorted_MDs:
			selection = temp[temp['MD'] == m] # get rows (read patterns) that have the given methylation density
			frac = np.nansum(selection['sample_C']) / float(total_C)
			weighted_fracs.append(frac)

		ctrl_weighted_MDs[ctrl] = weighted_fracs # append fractions to dataframe for given ctrl

	#--------------------------------------------------------------------------
	# reorder samples (columns) based on average methylation (meC of sample)
	reorder = zip(ctrlAveMeth, controls)
	reorder = sorted(reorder, key = lambda x: x[0])
	reorder_names = [i[1] for i in reorder]
	ctrl_weighted_MDs = ctrl_weighted_MDs[reorder_names]
	#--------------------------------------------------------------------------
	print 'Plotting weighted locus methylation density for controls...'
	# Plotting
	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])

	cmap = plt.get_cmap('coolwarm')
	colors = cmap(sorted_MDs) # colors based on the MDs, which should be between 0 and 1.0
	plot = ctrl_weighted_MDs.T.plot(kind='bar', stacked=True, ax=ax, color=colors, width=1, legend=True)

	plot.set_xlim( [-0.5, len(controls)-0.5] )
	ax.set_ylim(0,1)

	ax.set_yticks(np.arange(0,1.1,0.2))
	plt.yticks(fontsize=16)
	ax.set_ylabel('weighted sample fraction', fontsize=26)
	ax.yaxis.set_label_coords(-0.10, 0.5)

	# make bottom x-axis labels the number of reads
	ax.set_xticklabels([str('%.1e' % i) for i in ctrl_df[reorder_names].sum().values], rotation=90, fontsize=10) 
	ax.set_xlabel('total reads per sample input', fontsize=12)

	# colorbar methylation density legend
	# keep color scale from 0 to 1.0
	ax2 = fig.add_axes([0.1, 1.25, 0.5, 0.05])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1))
	cb.set_label('locus methylation density')
	cb.ax.xaxis.set_ticks_position('top')
	cb.set_clim(0.0, 1)

	# top x-axis labels = average methylation values
	if (vars(args)['averageMethylation'] == True):
		ax3 = ax.twiny()
		ax3.set_xticks(ax.get_xticks())
		ax3.set_xticklabels([round(i[0],3) for i in reorder], fontsize=10, rotation=45) 
		ax3.xaxis.set_ticks_position('top')
		ax3.xaxis.set_label_position('top')
		ax3.set_xlabel('average methylation', fontsize=12)
		ax3.xaxis.set_label_coords(0.5, 1.1)
		ax3.set_xlim(ax.get_xlim())

	# bottom label = sample type and number of samples
	ax4 = ax.twiny()
	ax4.set_xticks([])
	ax4.set_xticklabels([])
	ax4.set_xlabel(str(controlsType) + ' (n=' + str(len(controls)) + ')', fontsize=26)
	ax4.xaxis.set_label_coords(0.5, -0.25)

	ax.get_legend().remove()

	fig.savefig(fileTag + '.CTRLS_WEIGHTED_FRAC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)

	print 'Control Barplot created.'
	print 'Done'
	print ' '


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' readEFhistogram of sample set read fractions across samples '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

if (vars(args)['readEFhistogram'] == True):

	print "Generating Histogram of sample set read fractions across MD..."

	if inputFracAdjust == True:
		print 'Adjusted by relative sample input fraction'

	# get lists of MD values with occurances equal to number of reads with the given MD for each sample set

	# list for all reads from cases:
	cases_cols = ['MD'] + list(cases)
	# inputFrac Adjustment if needed:
	if inputFracAdjust == True:
		cases_md_sums = df[cases_cols].copy().set_index('MD').div(casesFracs).sum(axis=1) # sum of reads for each MD
	else:
		cases_md_sums = df[cases_cols].copy().set_index('MD').sum(axis=1) # sum of reads for each MD

	# list for all reads from controls:
	control_cols = ['MD'] + list(controls)
	# inputFrac Adjustment if needed:
	if inputFracAdjust == True:
		ctrl_md_sums = df[control_cols].copy().set_index('MD').div(ctrlFracs).sum(axis=1) # sum of reads for each MD
	else:
		ctrl_md_sums = df[control_cols].copy().set_index('MD').sum(axis=1) # sum of reads for each MD


	# adjustment for plotting in case the sums are actually fractions and thus the distr lists cannot be be made.
	# used to get list of whole number occurances of each MD
	# weights then adjusted to match this compensation such that the histogram EFs still match the actual read MD fractions across the samples of each sample set
	# essentially, the lowest non-zero fraction should be adjusted to be at least 10 or greater
	# as in, at least 10 more instances of the lowest nonzero fraction should be in the distr for plotting the histogram
	Adjust = int(round(10.0/np.min([i for i in list(cases_md_sums.values) + list(ctrl_md_sums.values) if i !=0.0])+1))


	cases_distr = [] # list to fill with MD occurances
	for loc, val in enumerate(cases_md_sums):
		cases_distr += int(val*Adjust) * [float(cases_md_sums.index[loc])]
	# list weights for histogram. Weights for all values is frequency of each read. ex: 100 reads --> 1/100 weight for each read
	cases_read_weights = (np.ones_like(cases_distr)/float(np.nansum(cases_md_sums)))/Adjust

	ctrl_distr = [] # list to fill with MD occurances
	for loc, val in enumerate(ctrl_md_sums):
		ctrl_distr += int(val*Adjust) * [float(ctrl_md_sums.index[loc])]
	# list weights for histogram. Weights for all values is frequency of each read. ex: 100 reads --> 1/100 weight for each read
	ctrl_read_weights = (np.ones_like(ctrl_distr)/float(np.nansum(ctrl_md_sums)))/Adjust


	# Sum of EFs of all MDs for a sample set should add up to 1.0

	#--------------------------------------------------------------------------

	# Plot Histogram

	f, ax = plt.subplots()

	# If # CpGs is less than 10, then would like each bin to represent
	# reads of a single MD. Except for the case of fully methylated reads.
	# These would end up in the last bin along with the second most methylated reads.
	# for bin edges, by default: [1,2), [2,3), [3,4), [4,5]
	# so for 0.125 MD increments (8 CpGs):
	# [0, 0.1249), [0.1249, 0.2499), [0.2499, 0.3749), [0.3749, 0.4999), [0.4999, 0.6249), [0.6249, 0.7499), [0.7499, 0.8749), [0.8749, 1.0]
	# 0 MD reads,  0.125 MD reads,   0.25 MD reads,    0.375 MD reads,   0.50 MD reads,    0.675 MD reads,   0.75 MD reads,    0.875 and 1.0 MD reads
	# Otherwise, bins can just span segments of 10% increments

	if max(df['numM'].values) >= 10:
		bins = 10
		xtick_pos = list(np.linspace(0.0, 1.0, bins, endpoint=False)) + [1.0]

	else:
		CpGs = max(df['numM'].values)
		xtick_pos = list(np.linspace(0.0, 1.0, CpGs, endpoint=False)) + [1.0]
		bins = [0.0] + [ i-0.0001 for i in xtick_pos[1:-1] ] + [1.0]

	cases_n, cases_bins, cases_patches = ax.hist(cases_distr, bins=bins, density=False, weights=cases_read_weights,
		color='red', align='mid', range=(0.0, 1.0), alpha=0.7, label=casesType + ' (n=' + str(len(cases)) + ')')

	ctrl_n, ctrl_bins, ctrl_patches = ax.hist(ctrl_distr, bins=bins, density=False, weights=ctrl_read_weights,
		color='skyblue', align='mid', range=(0.0, 1.0), alpha=0.6, label=controlsType + ' (n=' + str(len(controls)) + ')')

	# focus on EFs for methylated epialleles
	# background presumably mostly MD=0, so ignore and let go off axis
	heights = list(cases_n[1:]) + list(ctrl_n[1:]) # list of bar heights except the <10% MD background bars
	plt.ylim([0, max(heights)*1.25])
	ytick_pos = np.linspace(0.0, max(heights)*1.20, 5)

	# MD bins on x-axis:
	plt.xticks(xtick_pos, ['0%'] + [str(round(md*100,1))+'%' for md in xtick_pos[1:]], fontsize=18, rotation=30)

	# EF range on y-axis:
	plt.yticks(ytick_pos, ['0.0'] + ['%.1e' % i for i in ytick_pos[1:]], fontsize=18)

	ax.set_title('Fractions of reads', fontsize=22)
	plt.ylabel('sample set fraction', fontsize=20, labelpad=40)
	plt.xlabel('methylation density', fontsize=20, rotation=0)
	ax.yaxis.set_label_coords(-0.23,0.55)

	plt.legend(loc="upper center", fontsize=14, edgecolor='k', bbox_to_anchor=(0.55,0.90))
	plt.savefig(fileTag + '.ReadFrac-HIST.png', bbox_inches='tight', pad_inches=0.5, dpi=400)

	plt.close()

	print "Histogram completed"
	print ' '


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' TPR and FPR rates for Heat maps '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

print 'MDBC procedure:'
print ' '

fig = plt.figure() # for optional plotting within loop

for md in md_cutoffs:
	if md == 0.0:
		print 'ROC curve analysis for MD cutoff > ' + str(round(md, 2))
	else:
		print 'ROC curve analysis for MD cutoff >= ' + str(round(md, 2))

	# for each MD cutoff, compute weighted EF values for each sample.
	# weighted so that reads that cover only a few CpGs don't shift MD/EF significantly
	# whereas reads that cover many CpGs have larger effect because capture more of the methylation signature
	ctrl_weighted_ef_vals = []
	# compute EF for given MD cutoff for each control sample:
	for ctrl in controls:

		# total CpGs covered by all sample reads
		total_c_per_row = (df['numU'] + df['numM']) * df[ctrl]
		total_c = np.nansum(total_c_per_row)

		# select patterns with MDs at or above MD cutoff
		if md == 0:
			md_rows = df[df['MD'] > md] # if >= then EF would be 1.0 b/c all CpGs could be counted
		else:
			md_rows = df[df['MD'] >= md]

		# calculate number of CpGs covered by reads with MDs at or above cutoff for sample
		row_c = (md_rows['numU'] + md_rows['numM']) * md_rows[ctrl]
		sample_md_c = np.nansum(row_c)

		# determine EF or read counts of sample for given MD
		if usePeakCounts == True:
			sample_ef = (sample_md_c / numCpGs)
		else:
			sample_ef = (sample_md_c / float(total_c))

		# inputFrac adjustment if needed:
		if inputFracAdjust == True:
			fracAdjust = float(inputFracs[ctrl])
			weighted_ef = sample_ef / fracAdjust
		else:
			weighted_ef = sample_ef

		ctrl_weighted_ef_vals.append(weighted_ef)

	# for each MD cutoff, compute weighted EF values for each sample:
	case_weighted_ef_vals = []
	# compute EF for given MD cutoff for each case sample:
	for case in cases:

		# total CpGs covered by all sample reads
		total_c_per_row = (df['numU'] + df['numM']) * df[case]
		total_c = np.nansum(total_c_per_row)

		# select patterns with MDs at or above MD cutoff
		if md == 0:
			md_rows = df[df['MD'] > md]
		else:
			md_rows = df[df['MD'] >= md]

		# calculate number of CpGs covered by reads with MDs at or above cutoff for sample
		row_c = (md_rows['numU'] + md_rows['numM']) * md_rows[case]
		sample_md_c = np.nansum(row_c)

		# determine EF or read counts of sample for given MD
		if usePeakCounts == True:
			sample_ef = (sample_md_c / numCpGs)
		else:
			sample_ef = (sample_md_c / float(total_c))

		# inputFrac adjustment if needed:
		if inputFracAdjust == True:
			fracAdjust = float(inputFracs[case])
			weighted_ef = sample_ef / fracAdjust
		else:
			weighted_ef = sample_ef

		case_weighted_ef_vals.append(weighted_ef)


	# for each MD, generate ROC curves and find optimal EF cutoff
	ctrl_labels = [0 for i in ctrl_weighted_ef_vals]
	case_labels = [1 for i in case_weighted_ef_vals]
	roc_values = ctrl_weighted_ef_vals + case_weighted_ef_vals
	roc_labels = ctrl_labels + case_labels
	fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
	roc_auc = auc(fpr, tpr)

	# optimal EF threshold determined by maximum difference between TPR and FPR
	optimal_idx = np.argmax( tpr - fpr )
	optimal_threshold = thresholds[optimal_idx]
	sensitivity = tpr[optimal_idx]
	specificity = 1 - fpr[optimal_idx]
	
	# for each MD cutoff, record AUC, optimal EF cutoff, and associated sensitivity and specificity:
	aucs_per_md.append(roc_auc)
	opt_df_thresh_per_md.append(optimal_threshold)
	sens_per_md.append(sensitivity)
	spec_per_md.append(specificity)

	print "	AUC = " + str(round(roc_auc, 2))
	print "	optimal EF = " + str(optimal_threshold)
	print "	TPR = " + str(round(sensitivity, 2))
	print "	1-FPR = " + str(round(specificity, 2))


	# OPTIONAL: Generate ROC Curve for each MD cutoff:
	#--------------------------------------------------------------------------
	if (vars(args)['MDBCplots'] == True):
		print '	Generating ROC Curve...'

		lw = 2
		ax = fig.add_subplot(1,1,1)
		fig.set_size_inches(5, 4, forward=True)
		if md == 0.0:
			label = 'MDBC (MD > ' + str(round(md*100,2)) + '%' + '); AUC = %0.2f' % roc_auc
		else:
			label = 'MDBC (MD >= ' + str(round(md*100,2)) + '%' + '); AUC = %0.2f' % roc_auc
		plt.plot(fpr, tpr, color='red', lw=lw, label=label)

		plt.plot([0, 1], [0, 1], color='k', lw=lw, linestyle='--')
		plt.xlim([0.0, 1.00])
		plt.ylim([0.0, 1.00])
		plt.xlabel('FPR', fontsize=26)
		plt.ylabel('TPR', fontsize=26)
		#ax.yaxis.set_label_coords(-0.20,0.45)
		plt.xticks(np.arange(0,1.1,.2), [str(round(i,2)) for i in np.arange(0,1.1,.2)], fontsize=20)
		plt.yticks(np.arange(0,1.1,.2), [str(round(i,2)) for i in np.arange(0,1.1,.2)], fontsize=20)
		plt.legend(loc="lower right", fontsize=12, edgecolor='k')
		plt.savefig(fileTag + '.MD_' + str(round(md,2)) + '_cutoff-ROC.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
		plt.clf()
		print '	ROC Curve completed'
		print ' '

	a = case_weighted_ef_vals
	b = ctrl_weighted_ef_vals

	s, p = stats.ranksums(a, b)
	print '	# of cases = ' + str(len(a))
	print '	# of controls = ' + str(len(b))
	print '	Case MD cutoff median EF = ' + str(np.median(a))
	print '	Control MD cutoff median EF = ' + str(np.median(b))
	print '	ranksum p-val = ' + str(p)

	if (vars(args)['samples'] == True):
		print '	Generating Boxplot...'
		ax = fig.add_subplot(1,1,1)
		fig.set_size_inches(2, 4, forward=True)
		bp = ax.boxplot([a, b], positions=[1,3], widths=1, patch_artist=True)
		for flier in bp['fliers']:
			flier.set(marker='', color='black')
		for whisker in bp['whiskers']:
			whisker.set(color='black', linewidth=1)
		for cap in bp['caps']:
			cap.set(color='black', linewidth=1)
		for median in bp['medians']:
			median.set(color='black', linewidth=1)
		bp['boxes'][0].set( color='red', facecolor='salmon', linewidth=2, alpha=0.5)
		bp['boxes'][1].set( color='blue', facecolor='skyblue', linewidth=2, alpha=0.5)
		scatter = ax.scatter(x = np.random.normal(1, 0.1, size=len(a)),
			y = a, c='r', marker='.', edgecolors='', s = 50)
		scatter = ax.scatter(x = np.random.normal(3, 0.1, size=len(b)),
			y = b, c='b', marker='.', edgecolors='', s = 50)

		

		ax.axhline(y=optimal_threshold, linestyle='--', color='k')

		max_val = float(np.max(roc_values))
		if max_val != 0.0:
			s, p = stats.ranksums(a, b)
			if 0.01 <= p < 0.05:
				ax.text(x=1.7,y=max_val * 1.10, s='*', fontsize=30)
			if 0.001 <= p < 0.01:
				ax.text(x=1.6,y=max_val * 1.10, s='**', fontsize=30)
			if p < 0.001:
				ax.text(x=1.5,y=max_val * 1.10, s='***', fontsize=30)
			if p >= 0.05:
				ax.text(x=1.5,y=max_val * 1.12, s='ns', fontsize=30)
			ax.set_xlim([0, 4])
			ax.set_ylabel('epiallelic fraction', fontsize=20)
			#ax.yaxis.set_label_coords(-0.80,0.43)
			ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
			ax.set_ylim([0, max_val * 1.25])
			if md == 0.0:
				label = 'MD > ' + str(md)
			else:
				label = 'MD >= ' + str(md)
			ax.set_title(label, fontsize=20)
			ax.set_yticks(np.linspace(0, max_val * 1.15, 5))
			ax.set_yticklabels(['0.0'] + [ '%.1e' % i for i in np.linspace(0, max_val * 1.15, 5)[1:] ], fontsize=18)
		else:
			max_val = 0.0001
			s, p = stats.ranksums(a, b)
			if 0.01 <= p < 0.05:
				ax.text(x=1.7,y=max_val * 1.10, s='*', fontsize=30)
			if 0.001 <= p < 0.01:
				ax.text(x=1.6,y=max_val * 1.10, s='**', fontsize=30)
			if p < 0.001:
				ax.text(x=1.5,y=max_val * 1.10, s='***', fontsize=30)
			if p >= 0.05:
				ax.text(x=1.5,y=max_val * 1.12, s='ns', fontsize=30)
			ax.set_xlim([0, 4])
			ax.set_ylabel('epiallelic fraction', fontsize=20)
			#ax.yaxis.set_label_coords(-0.80,0.43)
			ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
			ax.set_ylim([0, max_val * 1.25])
			ax.set_title('MD = ' + str(md), fontsize=20)
			ax.set_yticks(np.linspace(0, max_val * 1.15, 5))
			ax.set_yticklabels(['0.0'] + [ '%.1e' % i for i in np.linspace(0, max_val * 1.15, 5)[1:] ], fontsize=18)

		fig.savefig(fileTag + '.MD_' + str(round(md,2)) + '_cutoff-BOXPLOT.png', bbox_inches='tight', pad_inches=1.0, dpi=400)
		plt.clf()

		print '	Boxplot completed'
	print ' '
	#--------------------------------------------------------------------------

	# Compute TPR and FPR for each EF and MD cutoff combination for TPR and FPR matrices

	# each row = a given MD cutoff, and each column (moving across row) = EF cutoff
	tprs = []
	fprs = []

	# Note that positives must be above the EF cutoff
	df_roc = pd.DataFrame({'label' : roc_labels, 'values' : roc_values})
	for ef in ef_cutoffs:
		true_positives = len(df_roc[ (df_roc['label'] == 1) & (df_roc['values'] > ef) ].index)
		false_negatives = len(df_roc[ (df_roc['label'] == 1) & (df_roc['values'] <= ef) ].index)
		true_negatives = len(df_roc[ (df_roc['label'] == 0) & (df_roc['values'] <= ef) ].index)
		false_positives = len(df_roc[ (df_roc['label'] == 0) & (df_roc['values'] > ef) ].index)
		tpr = float(true_positives)/(float(true_positives) + float(false_negatives))
		fpr = 1.0 - (float(true_negatives)/(float(true_negatives) + float(false_positives)))
		tprs.append(tpr)
		fprs.append(fpr)

	# append tpr/fpr rates of each EF to matrix row for the given md:
	tpr_matrix.append(tprs)
	fpr_matrix.append(fprs)

plt.close()

# save ROC information for each MD cutoff
roc_summary = pd.DataFrame()
roc_summary['MD'] = md_cutoffs
roc_summary['AUC'] = aucs_per_md
roc_summary['optimal_ef_cutoff'] = opt_df_thresh_per_md
roc_summary['TPR'] = sens_per_md
roc_summary['1-FPR'] = spec_per_md
roc_summary.to_csv(fileTag + '.SUMMARY.csv')


print "TPR and FPR matrices complete"
print ' '


#--------------------------------------------------------------------------

# Plotting Heatmaps

if (vars(args)['heatmaps'] == True):

	print "Generating Heat maps..."


	tpr_matrix = np.flipud(np.array(tpr_matrix))
	fpr_matrix = np.flipud(np.array(fpr_matrix))
	diff_matrix = tpr_matrix - fpr_matrix

	xtick_pos = np.arange(0, len(ef_cutoffs), 10)
	ytick_pos = np.arange(0, len(md_cutoffs), 4)



	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])
	cmap = plt.get_cmap('hot')
	ax.imshow(tpr_matrix, cmap=cmap, vmin=0, vmax=1) # option: remove vmin, vmax
	plt.xticks(xtick_pos, ['0.0'] + [ '%.1e' % ef for ef in ef_cutoffs[xtick_pos][1:] ], rotation=45, fontsize=14)
	plt.yticks(ytick_pos, [str(int(md*100))+'%' for md in md_cutoffs[::-1][ytick_pos]], fontsize=16)
	plt.xlabel('epiallelic fraction cutoff', fontsize=20)
	plt.ylabel('methylation\ndensity cutoff', fontsize=20, labelpad=50, rotation=0)
	ax.yaxis.set_label_coords(-0.26,0.3)
	ax2 = fig.add_axes([0.30, 0.77, 0.5, 0.03])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1)) #option: vmax = np.amax(tpr_matrix)
	cb.ax.set_title('TPR', fontsize=24)
	cb.ax.tick_params(labelsize=14) 
	cb.ax.xaxis.set_ticks_position('bottom')
	cb.set_clim(0.0, 1.0) # option: 0, np.amax(tpr_matrix)
	plt.savefig(fileTag + '.TPR-MATRIX.png', bbox_inches='tight', pad_inches=0.2, dpi=400)
	plt.close()



	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])
	cmap = plt.get_cmap('hot')
	ax.imshow(fpr_matrix, cmap=cmap, vmin=0, vmax=1) # option: remove vmin, vmax
	plt.xticks(xtick_pos, ['0.0'] + [ '%.1e' % ef for ef in ef_cutoffs[xtick_pos][1:] ], rotation=45, fontsize=14)
	plt.yticks(ytick_pos, [str(int(md*100))+'%' for md in md_cutoffs[::-1][ytick_pos]], fontsize=16)
	plt.xlabel('epiallelic fraction cutoff', fontsize=20)
	plt.ylabel('methylation\ndensity cutoff', fontsize=20, labelpad=50, rotation=0)
	ax.yaxis.set_label_coords(-0.26,0.3)
	ax2 = fig.add_axes([0.30, 0.77, 0.5, 0.03])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1)) #option: vmax = np.amax(tpr_matrix)
	cb.ax.set_title('FPR', fontsize=24)
	cb.ax.tick_params(labelsize=14) 
	cb.ax.xaxis.set_ticks_position('bottom')
	cb.set_clim(0.0, 1.0) # option: 0, np.amax(tpr_matrix)
	plt.savefig(fileTag + '.FPR-MATRIX.png', bbox_inches='tight', pad_inches=0.2, dpi=400)
	plt.close()



	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.9, 0.9])
	cmap = plt.get_cmap('hot')
	ax.imshow(diff_matrix, cmap=cmap, vmin=0, vmax=1) # option: remove vmin, vmax
	plt.xticks(xtick_pos, ['0.0'] + [ '%.1e' % ef for ef in ef_cutoffs[xtick_pos][1:] ], rotation=45, fontsize=14)
	plt.yticks(ytick_pos, [str(int(md*100))+'%' for md in md_cutoffs[::-1][ytick_pos]], fontsize=16)
	plt.xlabel('epiallelic fraction cutoff', fontsize=20)
	plt.ylabel('methylation\ndensity cutoff', fontsize=20, labelpad=50, rotation=0)
	ax.yaxis.set_label_coords(-0.26,0.3)
	ax2 = fig.add_axes([0.30, 0.77, 0.5, 0.03])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, orientation='horizontal', norm=mpl.colors.Normalize(vmin=0, vmax=1)) #option: vmax = np.amax(tpr_matrix)
	cb.ax.set_title('TPR - FPR', fontsize=24)
	cb.ax.tick_params(labelsize=14)
	cb.ax.xaxis.set_ticks_position('bottom')
	cb.set_clim(0.0, 1.0) # option: 0, np.amax(tpr_matrix)
	plt.savefig(fileTag + '.DIFF-MATRIX.png', bbox_inches='tight', pad_inches=0.2, dpi=400)
	plt.close()


	print "Heat maps completed"
	print ' '


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' ROC plots and Boxplots '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


# If roc_summary.csv already produced, comment the above code up to TPR and FPR rates and uncomment:
# roc_summary = pd.read_csv(sys.argv[4], index_col=0) # MDBCAnalysis_ZNF154reads_Lung_vs_Control_ROCsummaryTable.csv
# locus = str(sys.argv[4]).split('_')[1] # ex: MDBCAnalysis_ZNF154reads_Lung_vs_Control_ROCsummaryTable.csv --> ZNF154reads


# Overall optimal MD and EF defined as the largest positive difference of TPR - FPR for all MD/EF cutoffs tested:
diff_vals = roc_summary['TPR'].values - (1 - roc_summary['1-FPR'].values)
max_diff = np.max(diff_vals)
# in the case of ties, pick MD that also has highest AUC. If still a tie, then pick lowest MD cutoff.
max_idx = roc_summary.iloc[np.argwhere(diff_vals == max_diff).flatten()]['AUC'].idxmax()
opt_md = roc_summary.iloc[max_idx]['MD']
opt_ef = roc_summary.iloc[max_idx]['optimal_ef_cutoff']

print ' '
print 'Overall Optimal MD = ' + str(round(opt_md,2))
print 'Overall Optimal EF = ' + str(opt_ef)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' ROC curves for optimal MD/EF MDBC vs Average Methylation '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


# MDBC optimal MD sample EF values:

ctrl_weighted_ef_vals = []
# compute sample efs for optimal MD cutoff for each control sample:
for ctrl in controls:

	# total CpGs covered by all sample reads
	total_c_per_row = (df['numU'] + df['numM']) * df[ctrl]
	total_c = np.nansum(total_c_per_row)

	# select patterns with MDs at or above md cutoff
	if opt_md == 0:
		md_rows = df[df['MD'] > opt_md]
	else:
		md_rows = df[df['MD'] >= opt_md]

	# calculate number of CpGs covered by reads with MDs at or above cutoff for sample
	row_c = (md_rows['numU'] + md_rows['numM']) * md_rows[ctrl]
	sample_md_c = np.nansum(row_c)

	# determine EF or read counts of sample for given MD
	if usePeakCounts == True:
		sample_ef = (sample_md_c / numCpGs)
	else:
		sample_ef = (sample_md_c / float(total_c))

	# inputFrac adjustment if needed:
	if inputFracAdjust == True:
		fracAdjust = float(inputFracs[ctrl])
		weighted_ef = sample_ef / fracAdjust
	else:
		weighted_ef = sample_ef

	ctrl_weighted_ef_vals.append(weighted_ef)


case_weighted_ef_vals = []
# compute sample efs for optimal MD cutoff for each case sample:
for case in cases:

	# total CpGs covered by all sample reads
	total_c_per_row = (df['numU'] + df['numM']) * df[case]
	total_c = np.nansum(total_c_per_row)

	# select patterns with MDs at or above MD cutoff
	if opt_md == 0:
		md_rows = df[df['MD'] > opt_md]
	else:
		md_rows = df[df['MD'] >= opt_md]

	# calculate number of CpGs covered by reads with MDs at or above cutoff for sample
	row_c = (md_rows['numU'] + md_rows['numM']) * md_rows[case]
	sample_md_c = np.nansum(row_c)

	# determine EF or read counts of sample for given MD
	if usePeakCounts == True:
		sample_ef = (sample_md_c / numCpGs)
	else:
		sample_ef = (sample_md_c / float(total_c))

	# inputFrac adjustment if needed:
	if inputFracAdjust == True:
		fracAdjust = float(inputFracs[case])
		weighted_ef = sample_ef / fracAdjust
	else:
		weighted_ef = sample_ef

	case_weighted_ef_vals.append(weighted_ef)


# ROC case vs control Optimal MD EFs
ctrl_labels = [0 for i in ctrl_weighted_ef_vals]
case_labels = [1 for i in case_weighted_ef_vals]
roc_values = ctrl_weighted_ef_vals + case_weighted_ef_vals
roc_labels = ctrl_labels + case_labels
fpr_md, tpr_md, thresholds_md = roc_curve(roc_labels, roc_values)
roc_auc_md = auc(fpr_md, tpr_md)

# store values for CI estimation on ROC curve:
ci = pd.DataFrame({'labels':roc_labels, 'values':roc_values})
ci.to_csv(fileTag + '.OPTIMAL_MD_' + str(round(opt_md, 2)) + '_EF_VALS.csv')


#--------------------------------------------------------------------------

# Mean locus methylation

ctrl_ave_meth_values = []
for ctrl in controls:

	# total CpGs covered by all sample reads
	total_c_per_row = (df['numU'] + df['numM']) * df[ctrl]
	total_c = np.nansum(total_c_per_row)

	# total meCpGs covered by all sample reads
	total_meC_per_row = df['numM'] * df[ctrl]
	total_meC = np.nansum(total_meC_per_row)

	# compute average methylation (total number of meCpGs covered by fragments / total CpG's covered by fragments)
	ave_meC = total_meC / float(total_c)

	# inputFrac adjustment if needed:
	if inputFracAdjust == True:
		fracAdjust = float(inputFracs[ctrl])
		ave_meC = ave_meC / fracAdjust

	ctrl_ave_meth_values.append(ave_meC)


case_ave_meth_values = []
for case in cases:

	# total CpGs covered by all sample reads
	total_c_per_row = (df['numU'] + df['numM']) * df[case]
	total_c = np.nansum(total_c_per_row)

	# total meCpGs covered by all sample reads
	total_meC_per_row = df['numM'] * df[case]
	total_meC = np.nansum(total_meC_per_row)

	# compute average methylation (total number of meCpGs covered by fragments / total CpG's covered by fragments)
	ave_meC = total_meC / float(total_c)

	# inputFrac adjustment if needed:
	if inputFracAdjust == True:
		fracAdjust = float(inputFracs[case])
		ave_meC = ave_meC / fracAdjust

	case_ave_meth_values.append(ave_meC)


# ROC case vs control average methylation values
ctrl_labels = [0 for i in ctrl_ave_meth_values]
case_labels = [1 for i in case_ave_meth_values]
roc_values = ctrl_ave_meth_values + case_ave_meth_values
roc_labels = ctrl_labels + case_labels
fpr_ave, tpr_ave, thresholds_ave = roc_curve(roc_labels, roc_values)
roc_auc_ave = auc(fpr_ave, tpr_ave)

if (vars(args)['averageMethylation'] == True):
	# store values for CI estimation:
	ci = pd.DataFrame({'labels':roc_labels, 'values':roc_values})
	ci.to_csv(fileTag + '.aveMethValsForAUC-ci.csv')

#--------------------------------------------------------------------------

# ROC Curve:

fig, ax = plt.subplots()
lw = 2

# MDBC Optimal MD Cutoff ROC curve
print ' '
print 'Optimal MD cutoff ROC Curve:'
print '	MDBC (optimal MD cutoff = ' + str(round(opt_md,2)) + '):'
optimal_idx = np.argmax( tpr_md - fpr_md )
print '	AUC = ' + str(round(roc_auc_md, 2))
print '	optimal EF cutoff = ' + str(thresholds_md[optimal_idx])
print '	sensitivity = ' + str(tpr_md[optimal_idx])
print '	specificity = ' + str(1.0 - fpr_md[optimal_idx])

plt.plot(fpr_md, tpr_md, color='red',
         lw=lw, label='MDBC (MD = ' + str(round(opt_md*100,2)) + '%' + '); AUC = %0.2f' % roc_auc_md)


if (vars(args)['averageMethylation'] == True):
	# Average Methylation ROC curve
	print ' '
	print 'Versus Average Methylation ROC Curve:'
	if inputFracAdjust == True:
		print 'Adjusted by relative sample input fraction'
	optimal_idx = np.argmax( tpr_ave - fpr_ave )
	print '	AUC = ' + str(round(roc_auc_ave, 2))
	print '	optimal EF cutoff = ' + str(thresholds_ave[optimal_idx])
	print '	sensitivity = ' + str(tpr_ave[optimal_idx])
	print '	specificity = ' + str(1.0 - fpr_ave[optimal_idx])

	plt.plot(fpr_ave, tpr_ave, color='blue',
	         lw=lw, label='average methylation; AUC = %0.2f' % roc_auc_ave)


plt.plot([0, 1], [0, 1], color='k', lw=lw, linestyle='--')
plt.xlim([0.0, 1.00])
plt.ylim([0.0, 1.00])
plt.xlabel('FPR', fontsize=26)
plt.ylabel('TPR', rotation=0, fontsize=26, labelpad=40)
ax.yaxis.set_label_coords(-0.20,0.45)
plt.xticks(np.arange(0,1.1,.2), [str(round(i,2)) for i in np.arange(0,1.1,.2)], fontsize=20)
plt.yticks(np.arange(0,1.1,.2), [str(round(i,2)) for i in np.arange(0,1.1,.2)], fontsize=20)
plt.legend(loc="lower right", fontsize=12, edgecolor='k')

if (vars(args)['averageMethylation'] == True):
	rocTag = fileTag + '.OPTIMAL_MD_' + str(round(opt_md,2)) + '_vsAveMeth_ROC.png'
else:
	rocTag = fileTag + '.OPTIMAL_MD_' + str(round(opt_md,2)) + '_ROC.png'

plt.savefig(rocTag, bbox_inches='tight', pad_inches=0.5, dpi=400)
plt.close()


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' Box plots for optimal MD/EF MDBC vs Average Methylation '''

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# Optimal MD cutoff sample EF Boxplot

a = case_weighted_ef_vals
b = ctrl_weighted_ef_vals

fig, ax = plt.subplots(figsize=(2,4))
bp = ax.boxplot([a, b], positions=[1,3], widths=1, patch_artist=True)

for flier in bp['fliers']:
	flier.set(marker='', color='black')
for whisker in bp['whiskers']:
	whisker.set(color='black', linewidth=1)
for cap in bp['caps']:
	cap.set(color='black', linewidth=1)
for median in bp['medians']:
	median.set(color='black', linewidth=1)

bp['boxes'][0].set( color='red', facecolor='salmon', linewidth=2, alpha=0.5)
bp['boxes'][1].set( color='blue', facecolor='skyblue', linewidth=2, alpha=0.5)

scatter = ax.scatter(x = np.random.normal(1, 0.1, size=len(a)),
	y = a, c='r', marker='.', edgecolors='', s = 50)
scatter = ax.scatter(x = np.random.normal(3, 0.1, size=len(b)),
	y = b, c='b', marker='.', edgecolors='', s = 50)

s, p = stats.ranksums(a, b)
print ' '
print 'Boxplots Case vs Control using optimal MD cutoff:'
print '	# of cases = ' + str(len(a))
print '	# of controls = ' + str(len(b))
print '	Case median EF = ' + str(np.median(a))
print '	Control median EF = ' + str(np.median(b))
print '	ranksum p-val = ' + str(p)

# Identify optimal thresholds, plot, and print sens/spec:
case_labels = [1 for i in a]
ctrl_labels = [0 for i in b]
roc_values = a + b
roc_labels = case_labels + ctrl_labels
fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
roc_auc = auc(fpr, tpr)
optimal_idx = np.argmax( tpr - fpr )
optimal_threshold = thresholds[optimal_idx]
print '	optimal EF threshold = ' + str(optimal_threshold)


ax.axhline(y=optimal_threshold, linestyle='--', color='k')

max_val = np.max(roc_values)

s, p = stats.ranksums(a, b)
if 0.01 <= p < 0.05:
	plt.text(x=1.7,y=max_val * 1.10, s='*', fontsize=30)
if 0.001 <= p < 0.01:
	plt.text(x=1.6,y=max_val * 1.10, s='**', fontsize=30)
if p < 0.001:
	plt.text(x=1.5,y=max_val * 1.10, s='***', fontsize=30)
if p >= 0.05:
	plt.text(x=1.5,y=max_val * 1.12, s='ns', fontsize=30)


ax.set_xlim([0, 4])
ax.set_ylabel('epiallelic fraction', fontsize=18)
#ax.yaxis.set_label_coords(-0.97,0.43)
ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
ax.set_ylim([0, max_val * 1.25])
plt.title('MD >= ' + str(round(opt_md*100, 0)) + '%', fontsize=20)
plt.yticks(np.linspace(0, max_val * 1.15, 5), ['0.0'] + [ '%.1e' % i for i in np.linspace(0, max_val * 1.15, 5)[1:] ], fontsize=18)
fig.savefig(fileTag + '.OPTIMAL_MD_' + str(round(opt_md,2)) + '_BOXPLOT.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
plt.close()


#--------------------------------------------------------------------------

# Sample Average Methylation Boxplot

if (vars(args)['averageMethylation'] == True):

	a = case_ave_meth_values
	b = ctrl_ave_meth_values

	fig, ax = plt.subplots(figsize=(2,4))
	bp = ax.boxplot([a, b], positions=[1,3], widths=1, patch_artist=True)

	for flier in bp['fliers']:
		flier.set(marker='', color='black')
	for whisker in bp['whiskers']:
		whisker.set(color='black', linewidth=1)
	for cap in bp['caps']:
		cap.set(color='black', linewidth=1)
	for median in bp['medians']:
		median.set(color='black', linewidth=1)

	bp['boxes'][0].set( color='red', facecolor='salmon', linewidth=2, alpha=0.5)
	bp['boxes'][1].set( color='blue', facecolor='skyblue', linewidth=2, alpha=0.5)

	scatter = ax.scatter(x = np.random.normal(1, 0.1, size=len(a)),
		y = a, c='r', marker='.', edgecolors='', s = 50)
	scatter = ax.scatter(x = np.random.normal(3, 0.1, size=len(b)),
		y = b, c='b', marker='.', edgecolors='', s = 50)

	s, p = stats.ranksums(a, b)
	print ' '
	print 'Boxplots Case vs Control using optimal Average Methylation:'
	print '	# of cases = ' + str(len(a))
	print '	# of controls = ' + str(len(b))
	print '	Case median AveMeth = ' + str(np.median(a))
	print '	Control median AveMeth = ' + str(np.median(b))
	print '	ranksum p-val = ' + str(p)

	# Identify optimal thresholds, plot, and print sens/spec:
	case_labels = [1 for i in a]
	ctrl_labels = [0 for i in b]
	roc_values = a + b
	roc_labels = case_labels + ctrl_labels
	fpr, tpr, thresholds = roc_curve(roc_labels, roc_values)
	roc_auc = auc(fpr, tpr)
	optimal_idx = np.argmax( tpr - fpr )
	optimal_threshold = thresholds[optimal_idx]
	print '	optimal AveMeth threshold = ' + str(optimal_threshold)


	ax.axhline(y=optimal_threshold, linestyle='--', color='k')

	max_val = np.max(roc_values)

	s, p = stats.ranksums(a, b)
	if 0.01 <= p < 0.05:
		plt.text(x=1.7,y=0.90, s='*', fontsize=30)
	if 0.001 <= p < 0.01:
		plt.text(x=1.6,y=0.90, s='**', fontsize=30)
	if p < 0.001:
		plt.text(x=1.5,y=0.90, s='***', fontsize=30)
	if p >= 0.05:
		plt.text(x=1.5,y=0.95, s='ns', fontsize=30)

	ax.set_xlim([0, 4])
	ax.set_ylabel('average methylation', fontsize=18)
	#ax.yaxis.set_label_coords(-0.97,0.43)
	ax.set_xticklabels(['Cases', 'Controls'], fontsize=16, rotation=45)
	ax.set_ylim([0,1.1])
	plt.title('Average\nmethylation', fontsize=20)
	plt.yticks(np.arange(0,1.1,.2), [str(int(i*100)) + '%' for i in np.arange(0,1.1,.2)], fontsize=18)
	fig.savefig(fileTag + '.aveMeth-BOXPLOT.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
	plt.close()


print ' '
print 'MDBC ANALYSIS COMPLETED'
print ' '
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------






