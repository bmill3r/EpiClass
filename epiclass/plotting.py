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

import itertools
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.metrics import roc_curve, auc
from matplotlib import pyplot as plt
import matplotlib as mpl

#aesthetics for plots:
font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 12}
mpl.rc('font', **font)

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
        plt.xticks(self.xtick_pos, ['0%'] + [str(int(round(md*100, 0)))+'%' for md in self.xtick_pos[1:]],
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

        plt.xticks(self.xtickpos, ['0.0'] + ['%.1f' % x for x in self.xticks[self.xtickpos][1:]],
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
