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
import os
import os.path
from os import walk
import subprocess
from distutils.spawn import find_executable

import re
from collections import Counter
import numpy as np
import pandas as pd

from .logger import path_leaf

cwd = os.getcwd()

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
            print(' Check to make sure path is correct (Full Path Name) or samples have proper extensions.')
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
