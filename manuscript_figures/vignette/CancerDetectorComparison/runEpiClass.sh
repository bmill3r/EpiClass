#!/bin/bash

controlTrainSets="controlTrainingSets.txt"
caseTrainSets="caseTrainingSets.txt"

controlTestSets="controlTestSets.txt"
caseTestSets="caseTestSets.txt"

# make arrays to store sample set lists from each of the above files
# where each line is a set of samples for a given run
IFS=$'\n' read -d '' -r -a controlsTrain < $controlTrainSets
IFS=$'\n' read -d '' -r -a casesTrain < $caseTrainSets

IFS=$'\n' read -d '' -r -a controlsTest < $controlTestSets
IFS=$'\n' read -d '' -r -a casesTest < $caseTestSets

lociFile="_loci.txt"
IFS=$'\n' read -d '' -r -a loci < $lociFile
numLoci=$((${#loci[@]}-1))

fractions="sample_master_list.csv"


for i in {0..9} #10 runs, corresponds to number of Training/Test sets
do

	runNum=$(($i+1))

	echo "EpiClass Run "$runNum

	controlsSamplesTraining=${controlsTrain[i]}
	caseSamplesTraining=${casesTrain[i]}

	controlsSamplesTest=${controlsTest[i]}
	caseSamplesTest=${casesTest[i]}

	runTrainingDir="epiclass_MultiLoci_runs/Run_"$runNum"/training/"
	runTestingDir="epiclass_MultiLoci_runs/Run_"$runNum"/testing/"

	# for each run, perform epiclass on training and testing sets for each locus
	for l in $(seq 0 $numLoci)
	do

		locus=$(echo ${loci[l]} | cut -d' ' -f1)
		interval=$(echo ${loci[l]} | cut -d' ' -f2)
		locusREADtoDTpath=$(echo ${loci[l]} | cut -d' ' -f3) #READtoDT file generated based on a given locus genomic interval

		echo $locus
		
		locusTrainDirOut=$runTrainingDir$locus

		fileTagTrain=$locus"_Run_"$runNum"_training"
		locusTrainLogFile=$locusTrainDirOut"/"$fileTagTrain".log"

		mkdir -p $locusTrainDirOut

		echo "EpiClass training set..."

		epiclass MDBC -i $locusREADtoDTpath -a $caseSamplesTraining -b $controlsSamplesTraining \
		--fractions $fractions --ignoreEFsummary --optimalMDreadcounts \
		--fileTag $fileTagTrain -o $locusTrainDirOut > $locusTrainLogFile

		# get the optimal MD cutoff for read counts, obtained from the log file
		while read line; do
			if [[ "$line" == *"Optimal MD cutoff (read counts)"* ]]; then
	  			md=$(echo $line | cut -d'=' -f 2 | cut -d' ' -f 2)
			fi
		done < $locusTrainLogFile


		echo "MDBC test set with MD = "$md

		locusTestDirOut=$runTestingDir$locus

		fileTagTest=$locus"_Run_"$runNum"_testing"
		locusTestLogFile=$locusTestDirOut"/"$fileTagTest".log"

		mkdir -p $locusTestDirOut

		epiclass MDBC -i $locusREADtoDTpath -a $caseSamplesTest -b $controlsSamplesTest \
		--fractions $fractions --ignoreCountsummary --ignoreEFsummary --sampleValsAtMD $md \
		--fileTag $fileTagTest -o $locusTestDirOut > $locusTestLogFile

		# get training sample count values for optimal MD
		epiclass MDBC -i $locusREADtoDTpath -a $caseSamplesTraining -b $controlsSamplesTraining \
		--fractions $fractions --ignoreCountsummary --ignoreEFsummary --sampleValsAtMD $md \
		--fileTag $fileTagTrain -o $locusTrainDirOut

	done

done
