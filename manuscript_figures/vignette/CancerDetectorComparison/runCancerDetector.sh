#!/bin/bash

# run this script in the parent directory of "/inputReadFiles"
# for instance, in the directory: ./CancerDetectorComparisonZNF154/, which has inputReadFiles/,
# run: ../runCancerDetector.sh ZNF154

analysisLabel=$1 #ex: ZNF154, as the label, which will output files: reads.SAMPLE.run_#.tumor_burden
fileLabelLen=`echo $analysisLabel|wc -c`
fileLabelLen=$((${fileLabelLen} + ${1}))
echo $fileLabelLen

lambda=0.5 # a predefined lambda

mkdir -p cancerDetectorOutput

getLikelihoods="../CancerDetector/src/CalcReadLikelihood.py"
CancerDetector="../CancerDetector/src/CancerDetector.py"

readRunsDir=$(pwd)"/inputReadFiles/"

for run in $(ls "$readRunsDir")
do
	runNum=${run:4:2}
	echo "computing sample tumor burden for run "$runNum
	caseDir=$readRunsDir$run"/Cases/"
	controlDir=$readRunsDir$run"/Controls/"
	runMarkerFile=$readRunsDir$run"/"$analysisLabel".marker.run_"$runNum

	casesOutDir=cancerDetectorOutput/$run"/Cases/"
	controlsOutDir=cancerDetectorOutput/$run"/Controls/"

	mkdir -p $casesOutDir
	mkdir -p $controlsOutDir

	for control in $(ls "$controlDir")
	do
		echo $control
		name=${control:$fileLabelLen:26}
		readInput=$controlDir$control
		readLikelihood=$controlsOutDir$name".likelihood"
		readTumorBurden=$controlsOutDir$name".tumor_burden"

		python2.7 $getLikelihoods $readInput $runMarkerFile > $readLikelihood
		python2.7 $CancerDetector $readLikelihood $lambda > $readTumorBurden
	done

	for case in $(ls "$caseDir")
	do
		echo $case
		name=${case:$fileLabelLen:26}
		readInput=$caseDir$case
		readLikelihood=$casesOutDir$name".likelihood"
		readTumorBurden=$casesOutDir$name".tumor_burden"

		python2.7 $getLikelihoods $readInput $runMarkerFile > $readLikelihood
		python2.7 $CancerDetector $readLikelihood $lambda > $readTumorBurden
	done

done
