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

import numpy as np
import pandas as pd
import tables
from sklearn.metrics import roc_curve, auc

from .logger import path_leaf 

class mdbc():

    def __init__(self, df, cases, controls, fractions=None, mdcutoffs=None, maxCutoff=None, hdf_label='mdbc.mdvals.h5'):

        self.densityTable = pd.read_csv(df)
        self.mdvalues = sorted(list(set(self.densityTable['MD'].values)))

        self.maxCutoff = maxCutoff

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

        #if isinstance(fractions, pd.DataFrame):
        if fractions is not None:
            fractions = pd.read_csv(fractions)
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
                counts.append(np.nansum(selection))
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
        case_names = []
        for case in self.cases_filt:
            case_names.append(case)
            caseCounts.append(countCasesTable.loc[mdVals, case].sum())
            caseEFs.append(efCasesTable.loc[mdVals, case].sum())

        ctrlCounts = []
        ctrlEFs = []
        ctrl_names = []
        for ctrl in self.controls_filt:
            ctrl_names.append(ctrl)
            ctrlCounts.append(countControlsTable.loc[mdVals, ctrl].sum())
            ctrlEFs.append(efControlsTable.loc[mdVals, ctrl].sum())

        countVals = pd.DataFrame(
            {'cases': pd.Series(caseCounts), 'controls': pd.Series(ctrlCounts),
            'case_samples': pd.Series(case_names), 'control_samples': pd.Series(ctrl_names)})
        EFVals = pd.DataFrame(
            {'cases': pd.Series(caseEFs), 'controls': pd.Series(ctrlEFs),
            'case_samples': pd.Series(case_names), 'control_samples': pd.Series(ctrl_names)})

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

        if self.maxCutoff is None:
            # use non-MD 0% values
            max_ef = pd.concat(self.readCountsPerMDtables, axis=1).loc[self.mdvalues[1:]].max().max() * 1.10
        else:
            max_ef = self.maxCutoff
            #max_ef = 100
        
        return np.array(list(np.linspace(0.000, max_ef, 100, endpoint=False)) + [max_ef])

    @property
    def readEFCutoffRange(self):
        '''
        Define the range of possible read cutoffs to use for the fraction of methylated reads in a sample
        to classiy it as positive, given a MD cutoff.
        '''

        if self.maxCutoff is None:
            # use non-MD 0% values
            max_ef = pd.concat(self.readEFsPerMDtables, axis=1).loc[self.mdvalues[1:]].max().max() * 1.10
        else:
            max_ef = self.maxCutoff

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
