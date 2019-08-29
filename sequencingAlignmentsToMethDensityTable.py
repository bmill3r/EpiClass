#!/usr/bin/env	python

import pandas as pd
import numpy as np
import sys
import os
from os import walk
import subprocess
from distutils.spawn import find_executable
import re
from collections import Counter
import datetime
import ntpath
import argparse
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
	description=
	'''
	Convert reads in sequencing alignment files (SAM/BAM) from genomic intervals of interest into a methylationDensityTable.csv.

	Alignment files are sorted alignment files from Bismark. If BAM, requires samtools view to access the reads.
	Otherwiise, will need alignment files to be in SAM.

	Extracted Methylation sequences come from reads that overlap intervals in interest. Sequence can be entire sequence or only sequence in overlap from a read or stitched together from paired-reads.

	example of a line (read) from a bismark alignment file that is expected:
	readID 163 chrom readStart 40 120M = 74960987 158 actualDNASequence alignmentScore NM:i:30 MD:Z:5G8G0G5G0G16G0G4G1G3G1G6G0G3G6G2G0G0T0G3G2G0G4G4G0G1G5G3G1G6G1 XM:Z:methylationInfo(...z...h...x..Z...etc) XR:Z:GA XG:Z:GA XB:Z:6

	Methylation information taken from the 'XM:Z:' line position.

	Requires a directory that contains the alignment files of interest. Each SAM/BAM file will be considered a separate sample to include in the methylationDensityTable.csv.

	Requires a tab-delimited BED file where each row contains the chromosome, start, and stop coordinates of interval(s) of interest to extract reads from.
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

	In addition to the methylationDensityTable.csv, returns a record .out file with STDOUT text.
	''',

	epilog=
	'''
	usage:
	./sequencingAlignmentsToMethDensityTable.py -in pathToFiles/ -i intervals.bed

	./sequencingAlignmentsToMethDensityTable.py -in inputDir/ -i chr19:58220244-58220744 -rs -out /Users/home/project/ -ov 0.75 -t ZGTreads10kb
	(only keep reads that have at least 75% overlap with the specific interval and only extract read methylation sequences that are contained within the overlap.)
	(also outputs files to /Users/home/project/, and only consideres alignment files in inputDir/ that contain 'ZGTreads10kb' in their name. )
	'''
	)


parser.add_argument('-in', '--inputDir', metavar='pathToFiles/', required=True, action="store",
	help='''REQUIRED.
	Path to directory containing the files of interest.
	Each sample of interest should have its own read file (.sam, for instance) in the directory.
	Path should end with a '/' (for Mac/Unix)
	'''
	)

parser.add_argument('-t', '--fileType', metavar='fileType', action="store",
	help='''
	Label indicating the type of files to select.
	Typically contains information like the locus.
	For example: 'ZNF154reads.sam' will indicate selection of the *ZNF154reads.sam* files in the input directory.
	i.g., HCC1.dedup.SOcoord.ZNF154reads.sam (sample = HCC1)
	'''
	)

parser.add_argument('-i', '--interval', metavar='chr:start-stop', required=True,
	help='''REQUIRED.
	Coordinates of interval of interest to select reads from.
	Can also be a tab-delimited 'intervalFile.txt' or BED file with no column names.
	Coordinates of intervals organized in each row as follows:
	chr    start    stop
	'''
	)

parser.add_argument('-ov', '--overlap', metavar='overlap', action="store", type=float,
	help='''
	Extract methylation sequences from read/read pairs that overlap the defined interval at least this much.
	'''
	)

parser.add_argument('-rs', '--readSlice', action='store_true', default=False,
	help='''
	Flag call to use methylation sequences information only from sections of selected reads that are contained in the overlap with the interval(s).
	Otherwise use methylation information from the entire read(s) that overlap.
	'''
	)

parser.add_argument('-out', '--outDir', action='store',
	help='''
	Directory path to store output files.
	'''
	)

#----------------------------------------------------------------------------------------

# get file name from path no matter the operating system or file path
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#----------------------------------------------------------------------------------------
# Set up variables:

# Collect command line arguments
args = parser.parse_args()

# Date run:
date = datetime.date.today()

# current working directory
cwd = os.getcwd()

# make sure input directory ends with back slash:
if str(args.inputDir).endswith('/'):
	path = args.inputDir
else:
	path = str(args.inputDir) + '/'


# check if interval file used, otherwise use single interval input of chr:start-stop
if os.path.exists(args.interval) == True:
	intervalTag = path_leaf(args.interval)
	# should be tab-delimited with no column names.
	# Essentially a BED file with intervals organized as chr, start, stop
	intervalFileName = args.interval
	intervals = pd.read_table(intervalFileName, names=['chr', 'start', 'stop'], index_col=False)
else:
	# chr:start-stop
	intervalTag = args.interval
	int_chr = args.interval.split(':')[0]
	int_start = int(args.interval.split(':')[1].split('-')[0])
	int_stop = int(args.interval.split(':')[1].split('-')[1])
	intervals = pd.DataFrame({'chr':[int_chr], 'start':[int_start], 'stop':[int_stop]})
	# save as temporary bedFile
	intervalFileName = str(cwd) + '/' + 'interval.txt'
	intervals.to_csv(intervalFileName, columns = ['chr', 'start', 'stop'], header=False, index=False, sep='\t')
	intervalTag = args.interval


if vars(args)['overlap'] != None:
	overlapTag = args.overlap
else:
	overlapTag = 'ANY'


# tags to name files, incorporating output directory if necessary
if vars(args)['outDir'] != None:
	outdir = str(args.outDir)
	fileTag = outdir + 'sequencingAlignmentsToMethDensityTable.' + str(args.fileType) + '.Interval=' + intervalTag + '.OverlapThresh=' + str(overlapTag) + '.readSlice=' + str(args.readSlice) + '.Created=' + str(date)
else:
	fileTag = 'sequencingAlignmentsToMethDensityTable.' + str(args.fileType) + '.Interval=' + intervalTag + '.OverlapThresh=' + str(overlapTag) + '.readSlice=' + str(args.readSlice) + '.Created=' + str(date)


# print std.out 'print' statements to new file AND to the console:
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


#----------------------------------------------------------------------------------------

print 'Date:' + str(date)
print '#--------------------------------------------------------------------------'
print 'Command:'
print sys.argv
print ' '
print args
print '#--------------------------------------------------------------------------'
print 'Building Methylation Density Table for samples of fileType: ' + str(args.fileType) + '...'
print 'Locating files in ' + path
print 'Selecting reads from each file with ' + str(overlapTag) + ' overlap with interval(s):'
print str(args.interval)

#----------------------------------------------------------------------------------------

# final dataframe to populate
# contains read counts per each sample for all reads which overlap any of the intervals
# also include chr, and start/end coordinates of the selected read methylation information.
merged_df = pd.DataFrame()


files = []
for (dirpath, dirnames, filenames) in walk(path):
    files.extend(filenames)

    # filer for .bam/.sam files:
    files = [i for i in files if bool(re.search(r'[sb]am$', i, re.IGNORECASE)) ]

    # only process files of the given type
    if vars(args)['fileType'] != None:
		files = [i for i in files if str(args.fileType) in i]

    for f in files:

		# append reads and corresponding methylation information for reads that pass cutoffs for each sample
		sample_reads = []

		sample_name = f.split('.', 1)[0] # expects sample same to be the first string before sep in file name.
		inputFile = path + f
		print ' '
		print 'extracting reads from: ' + path_leaf(f)

		lines = [] # collect reads from file
		#----------------------------------------------------------------------------------------
		# check if samtools exists:
		if find_executable('samtools') is not None:

			# if samtools exists, use to read file and select reads overlapping intervals:
			cmd = 'samtools view -L ' + str(intervalFileName) + ' ' + str(inputFile)
			print ' accessing reads with `samtools view`.'
			print ' ' + cmd

			# run code, check for errors, get output
			reads = subprocess.Popen([cmd], shell=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)
			stdout, stderr = reads.communicate()
			# if error, end script.
			if len(stderr) != 0:
				print 'Error with `samtools view`.'
				print stderr
				print 'Exiting script.'
				sys.exit(0)
			# otherwise collect lines (reads) of interest:
			else:
				for line in stdout.splitlines():
					lines.append([ line.split('\t')[i] for i in [0,2,3,13] ])


		# if samtools does not exist, and sample is sam:
		elif bool(re.search(r'sam$', inputFile, re.IGNORECASE)):

			with open(inputFile, 'r') as sam:
				mylist = sam.read().splitlines()
				for line in mylist:
					# ignore header lines
					if line.startswith('@'): 
						continue
					else:
						# read_name, chr, read_start, meth
						lines.append([ line.split('\t')[i] for i in [0,2,3,13] ])


		# if samtools does not exist, and sample is bam:
		elif bool(re.search(r'bam$', inputFile, re.IGNORECASE)):
			print ' `samtools view` needed to access reads in BAM file.'
			print ' skipping file...'
			continue

		#----------------------------------------------------------------------------------------
		df = pd.DataFrame(lines, columns=['read', 'chr', 'read_start', 'meth'])
		df['meth'] = df.copy()['meth'].str.split(':', 2).str[2] # remove 'XM:Z:' at beginning of meth sequence
		df['readLength'] = df['meth'].str.len()
		df['read_start'] = df['read_start'].astype(int)
		df['read_end'] = df['read_start'] + df['readLength']

		# get selected reads for each interval (doesn't matter if paired or not)
		# any paired reads in which both pairs overlap, methylation information will be stitched together.
		# otherwise, methylation information for the individual reads will be extracted.
		for interval in intervals.itertuples():

			interval_reads = [] # collect methylation information for sample reads that overlap interval

			int_chr = interval.chr
			int_start = interval.start
			int_stop = interval.stop
			intervalname = str(interval.chr) + ':' + str(int(interval.start)) + '-' + str(int(interval.stop))

			# keep reads that are on same chromosome and overlap interval
			readsIn = df[ (df['chr'] == int_chr) ]
			readsIn = readsIn[ ~(int_start >= readsIn['read_end']) ]
			readsIn = readsIn[ ~(readsIn['read_start'] >= int_stop) ]
			readsIn['int_start'] = int_start
			readsIn['int_stop'] = int_stop

			# get length and proportion of read overlap to corresponding interal

			Ecols = readsIn[['read_end', 'int_stop']].idxmin(axis=1) # overlap end
			readsIn['E'] = [ readsIn.iloc[i][Ecols.iloc[i]] for i in np.arange(0,len(Ecols)) ]

			Scols = readsIn[['read_start', 'int_start']].idxmax(axis=1) # overlap start
			readsIn['S'] = [ readsIn.iloc[i][Scols.iloc[i]] for i in np.arange(0,len(Scols)) ]

			readsIn['overlap'] = (readsIn['E'] - readsIn['S']) / readsIn['readLength']
			readsIn['sliceLength'] = (readsIn['overlap'] * readsIn['readLength']).astype(int)
			readsIn['outsideLength'] = readsIn['readLength'] - readsIn['sliceLength']

			# if --overlap, keep reads that overlap the interval by at least this fraction
			if vars(args)['overlap'] != None:
				readsIn = readsIn[readsIn['overlap'] >= args.overlap]

			print ' ' + str(len(readsIn)) + ' reads/read-pairs overlapping ' + intervalname

			# choose whether to use methylation information that is only contained in the read overlap, or include the entire read
			if vars(args)['readSlice'] == False:
				s = 'meth'
			else:
				s = 'slice'
				# get slices of meth seq contained in overlap:
				slices = []
				for index, row in readsIn.iterrows():
					sliceStart = row['S'] - row['read_start']
					sliceEnd = len(row['meth']) - (row['read_end'] - row['E'])
					slices.append(row['meth'][sliceStart:sliceEnd])
				readsIn[s] = slices
			
			# combine paired-reads and stitch together methylation pattern and extract
			for g, d in readsIn.groupby('read'):

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
							stitchedMeth = d[s].iloc[0] + d[s].iloc[1][readOverlap:]
						if s == 'meth':
							readOverlap = d['read_end'].iloc[0] - d['read_start'].iloc[1]
							stitchedMeth = d[s].iloc[0] + d[s].iloc[1][readOverlap:]

					# if end of R2 overlaps beginning of R1:
					if (d['read_end'].iloc[1] > d['read_start'].iloc[0]) and (d['read_start'].iloc[0] > d['read_start'].iloc[1]):
						if s == 'slice':
							readOverlap = d['E'].iloc[1] - d['S'].iloc[0]
							stitchedMeth = d[s].iloc[1] + d[s].iloc[0][readOverlap:]
						if s == 'meth':
							readOverlap = d['read_end'].iloc[1] - d['read_start'].iloc[0]
							stitchedMeth = d[s].iloc[1] + d[s].iloc[0][readOverlap:]

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
							methSeqStart = np.min([ d['S'].iloc[0], d['S'].iloc[1] ])
							methSeqStop = np.max([ d['E'].iloc[0], d['E'].iloc[1] ])
						if s == 'meth':
							methSeqStart = np.min([ d['read_start'].iloc[0], d['read_start'].iloc[1] ])
							methSeqStop = np.max([ d['read_end'].iloc[0], d['read_end'].iloc[1] ])
					else:
						if s == 'slice':
							methSeqStart = d['S'].iloc[0]
							methSeqStop = d['E'].iloc[0]
						if s == 'meth':
							methSeqStart = d['read_start'].iloc[0]
							methSeqStop = d['read_end'].iloc[0]

					interval_reads.append([chrom, d['int_start'].iloc[0], d['int_stop'].iloc[0], methSeqStart, methSeqStop, len(stitchedMeth), numU, numM, md])

			if len(interval_reads) == 0:
				print '	Warning: No reads found overlapping with ' + intervalname + ' in ' + str(inputFile)
				continue
			else:
				print ' ' + str(len(interval_reads)) + ' unique sequences with CpGs'
				for i in interval_reads:
					sample_reads.append(i)

		# check to see if any reads were selected:
		if len(sample_reads) == 0:
			print '	Warning: No reads found overlapping any intervals in: ' + str(inputFile)
			continue
		else:
			print ' ' + str(len(sample_reads)) + ' total unique sequences with CpGs'
			# Build dataframe from sample_reads array for selected sample reads only
			cols = ['chr', 'interval_start', 'interval_stop', 'methSeqStart', 'methSeqStop', 'methSeqLength', 'numU', 'numM', 'MD']
			sample_reads_df = pd.DataFrame(sample_reads, columns=cols)
			sample_reads_df.set_index(cols, inplace=True)

			# get counts of unique reads
			unique_read_counts = dict(Counter(sample_reads_df.index.tolist())) 

			# index becomes the read, and the column becomes the number of that particular read in the sample
			unique_read_counts_df = pd.DataFrame.from_dict(unique_read_counts, orient='index')

			# if same sample name had been used, append a counter to the end
			if sample_name in merged_df.columns:
				sample_name = sample_name + '.' + str( len( [i for i in merged_df.columns.tolist() if sample_name in i] ) )

			# remake dataframe with counts
			unique_read_counts_df.columns = [sample_name]

			# append the unique_read_counts_df sample dataframe to the final dataframe
			merged_df = merged_df.join(unique_read_counts_df, how='outer')

print ' '
print 'Cleaning table...'
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
merged_df.sort_values(['chr', 'methSeqStart', 'methSeqStop'], inplace=True)

# save final dataframe
merged_df.to_csv(fileTag + '.csv', index=False )

# remove temporary interval file, if exists:
if os.path.exists(str(cwd) + '/' + 'interval.txt') == True:
	os.remove(str(cwd) + '/' + 'interval.txt')

print 'Finished'
print ' '


