#!/usr/bin/env	python

import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as Colors
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
import scipy
import sklearn
from sklearn.metrics import roc_curve, auc, roc_auc_score
import scipy.stats as stats


font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 12}
matplotlib.rc('font', **font)

#------------------------------------------------------------------------------------
'''
./Figure1_average_methylation_lineplot.py
	Figure1_OV_FT_Uterine_ZNF154_probes.csv
	Figure1_WBC_ZNF154_probes.csv
	Figure1_TCGA_450K_ZNF154_probes.csv
'''



probes = ['cg11294513', 'cg05661282', 'cg21790626', 'cg27049766']
x_pos = [0.0, 1.35, 3.58, 4.0] # roughly based on relative bp distance between probes


d1 = pd.read_csv(sys.argv[1]) # Figure1_OV_FT_Uterine_ZNF154_probes.csv --> ft, ucec, ov probes

fts = d1[[col for col in d1.columns if 'ft' in col]]
fts['cg11294513'] = fts[[col for col in fts.columns if 'cg11294513' in col]]
fts['cg05661282'] = fts[[col for col in fts.columns if 'cg05661282' in col]]
fts['cg21790626'] = fts[[col for col in fts.columns if 'cg21790626' in col]]
fts['cg27049766'] = fts[[col for col in fts.columns if 'cg27049766' in col]]
fts = fts[[col for col in fts.columns if 'ft' not in col]]

ucecs = d1[[col for col in d1.columns if 'ucec' in col]]
ucecs['cg11294513'] = ucecs[[col for col in ucecs.columns if 'cg11294513' in col]]
ucecs['cg05661282'] = ucecs[[col for col in ucecs.columns if 'cg05661282' in col]]
ucecs['cg21790626'] = ucecs[[col for col in ucecs.columns if 'cg21790626' in col]]
ucecs['cg27049766'] = ucecs[[col for col in ucecs.columns if 'cg27049766' in col]]
ucecs = ucecs[[col for col in ucecs.columns if 'ucec' not in col]]

ovs = d1[[col for col in d1.columns if 'ov' in col]]
ovs['cg11294513'] = ovs[[col for col in ovs.columns if 'cg11294513' in col]]
ovs['cg05661282'] = ovs[[col for col in ovs.columns if 'cg05661282' in col]]
ovs['cg21790626'] = ovs[[col for col in ovs.columns if 'cg21790626' in col]]
ovs['cg27049766'] = ovs[[col for col in ovs.columns if 'cg27049766' in col]]
ovs = ovs[[col for col in ovs.columns if 'ov' not in col]]

wbcs = pd.read_csv(sys.argv[2], index_col=0) # Figure1_WBC_ZNF154_probes.csv


ft_aves = fts.mean(axis=0)
ucec_aves = ucecs.mean(axis=0)
ov_aves = ovs.mean(axis=0)
wbc_aves = wbcs.mean(axis=1)


#------------------------------------------------------------------------------------
''' plot averages, lineplot '''


fig, ax = plt.subplots(figsize=(6.5,2))

#ax.plot(x_pos, ov_aves.values, color='red', label='Ovarian carcinoma')
#ax.plot(x_pos, ucec_aves.values, color='navy', label='Uterine')
#ax.plot(x_pos, ft_aves.values, color='blue', label='Fallopian tube')
#ax.plot(x_pos, wbc_aves.values, color='deepskyblue', label='WBC')

# cycle colors:
color_blues=iter(cmx.winter(np.linspace(0,1,12)))
color_reds=iter(cmx.autumn(np.linspace(0,1,14)))


tissues = ['BLCA', 'COAD', 'HNSC', 'LIHC', 'LUSC', 'PRAD', 'STAD', 'UCEC']
df2 = pd.read_csv(sys.argv[3], index_col=0) # Figure1_TCGA_450K_ZNF154_probes.csv


# plot tumors
ct=next(color_reds)
ax.plot(x_pos, ov_aves.values, color=ct, label='EOCs *')

for tissue in tissues:
	temp_df = df2[[i for i in df2.columns if tissue in i]]
	tumors = temp_df[[i for i in temp_df.columns if '_tumors' in i]].values[:4]
	ct=next(color_reds)
	ax.plot(x_pos, tumors, c=ct, label=tissue)
	

# plot controls
cb=next(color_blues)
ax.plot(x_pos, wbc_aves.values, color=cb, label='WBCs **')

for tissue in tissues:
	temp_df = df2[[i for i in df2.columns if tissue in i]]
	ctrls = temp_df[[i for i in temp_df.columns if '_ctrls' in i]].values[:4]
	cb=next(color_blues)
	ax.plot(x_pos, ctrls, c=cb, label=tissue + ' ')



# aesthetics and labels
plt.xticks(x_pos, probes, fontsize=14, rotation=45)
plt.yticks(np.arange(0,1.1,0.2), [str(int(i*100)) + '%' for i in np.arange(0,1.1,0.2)], fontsize=14)

ax.set_ylabel('average\nmethylation', rotation=0, labelpad=40, fontsize=20)
ax.yaxis.set_label_coords(-0.31,0.3)

plt.legend(bbox_to_anchor=(1.04,1.15), loc="upper left", fontsize=10, edgecolor='k', ncol=2)

plt.text(x=4.5, y=1.15, s='Tumors')
plt.text(x=5.5, y=1.15, s='Controls')

fig.savefig('Figure1_average_methylation_lineplot.png', bbox_inches='tight', pad_inches=0.5, dpi=600)




















