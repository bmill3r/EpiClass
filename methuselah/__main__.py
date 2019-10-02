'''
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

'''

"""IMPORT DEPENDENCIES"""

import sys
import os
import datetime

import pandas as pd
from matplotlib import pyplot as plt

"""IMPORT MODULES"""
from .arguments import get_arguments
from .logging import Logger, path_leaf
from .reading import dreamingToDensityTable, readsToDensityTable
from .analyzing import mdbc
from .plotting import boxplot, boxplot2sets, stackedBarplot, histogram, heatmap, rocplot

def main(args=None):

	now = datetime.datetime.now()
	timestamp = now.strftime("%Y-%m-%d_%H-%M")
	date = str(datetime.date.today())
	cwd = os.getcwd()

	sys.stdout = Logger()

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



