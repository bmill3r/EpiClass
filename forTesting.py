#!/usr/bin/env  python

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

Optimizing and predicting performance of DNA methylation biomarkers using methylation density information.

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

# /Users/millerbf2/Desktop/Elnitski/Projects/Methylation_Density_Genomewide/methuselah
# ./forTesting.py MDBC -i ../MDBC_Analysis/OV_plasma/DREAMtoMD.DT.20190211-DREAMing_well_melt_temps_raw.input2bg.csv -a JKLR FLAB JGLR JFGH GKHG JHPP WHPA HKLP FQTR JUTM TFLK JQRT WQHT QWID QWER POIU ASDF GFUR ZXCV HJKL REWQ SDFG XCVB JGIS WWZX KRUK -b 101086 102801 100626 101425 145355 100296 103197 100250 102598 101599 100292 106732 106853 109195 109286 109336 110164 113493 114250 121203 121274 123154 123874 124110 127216 128141 129125 129499 130752 131004 139606 103971 101997 102073 106401 106136 109837 121624 102060 103782 114482 109604 101968 102919 100654


# IMPORT DEPENDENCIES
import sys
import os
import datetime
import pandas as pd
from matplotlib import pyplot as plt
# IMPORT MODULES
from EpiClass.arguments import get_arguments
from EpiClass.logger import Logger, path_leaf, Logging
from EpiClass.reading import dreamingToDensityTable, readsToDensityTable
from EpiClass.analyzing import mdbc
from EpiClass.plotting import boxplot, boxplot2sets, stackedBarplot, histogram, heatmap, rocplot

now = datetime.datetime.now()
timestamp = now.strftime("%Y-%m-%d_%H-%M")
date = str(datetime.date.today())
cwd = os.getcwd()
args, args_dict = get_arguments(args=None)

if args_dict['cmd'] == 'MDBC':
    df = args_dict['input']
    cases = args_dict['cases']
    controls = args_dict['controls']
    fractions = args_dict['fractions']
    mdcutoffs = args_dict['MDcutoffs']
    outdir = args_dict['output']
    hdf_label = args_dict['hdf_label']

    dfName = path_leaf(df).split('.DT.')[1].split('.csv')[0]

    if outdir is not None:
        if outdir.endswith('/'):
            pass
        else:
            outdir = outdir + '/'
        filename = outdir + 'MDBC.' + dfName
    else:
        filename = cwd + '/MDBC.' + dfName

    if hdf_label is not None:
        pass
    else:
        hdf_label = cwd + '/MDBC_H5DF.' + dfName + '.h5'

    if args_dict['fileTag'] is not None:
        fileTag = args_dict['fileTag']
        filename = filename + '.' + fileTag


    classifier = mdbc(df=df, cases=cases, controls=controls,
                     fractions=fractions, mdcutoffs=mdcutoffs,
                     hdf_label=hdf_label)

    print(classifier.inputFracs)


#print(classifier.readEFCutoffRange)

#print(pd.concat([classifier.CpGcolumns, classifier.normalizedDensTable], axis=1))

#print(classifier.readEFsPerMDtables)

#print(pd.concat(classifier.readEFsPerMDtables, axis=1).loc[classifier.mdvalues[1:]])








