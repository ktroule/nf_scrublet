#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -- REQUIRED LIBRARIES
import matplotlib.pyplot as plt
import scrublet as scr
import seaborn as sns
import scanpy as sc
import numpy as np
import sys, getopt
import argparse
import scipy
import os
import re

# -- CUSTOM MODULES
import modules.pval_correction as kcorr

# -- DEFINE SCRIPT INPUT ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample')
parser.add_argument('-d', '--directory')
parser.add_argument('-g', '--min_genes')
parser.add_argument('-c', '--min_cells')
parser.add_argument('-r', '--rnd_seed')
parser.add_argument('-out', '--out_dir')

args = parser.parse_args()
print(args)


# -- SCRUBLET SCORE CALCULATION
adata = sc.read_10x_mtx(os.path.join(args.directory,
                                     args.sample,
                                     'filtered_feature_bc_matrix'),
                        cache = True)

adata.var_names_make_unique()

# -- Rename Sample barcode by adding sample id
adata.obs_names = [ args.sample + '_' + i for i in adata.obs_names ]

# -- Basic initial filtering
sc.pp.filter_cells(adata,
                    min_genes = int(args.min_genes))

sc.pp.filter_genes(adata,
                   min_cells = int(args.min_cells))

sc.pp.calculate_qc_metrics(adata,
                            inplace = True,
                            percent_top = None)

# -- Calculate doublet score with scrublet
np.random.seed(int(args.rnd_seed))

scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose = False)
adata.obs['scrublet_score'] = doublet_scores
adata.obs['scrublet_prediction'] = predicted_doublets

# -- Plot scrublet distribution
sns.displot(adata.obs['scrublet_score']).set(title = args.sample)
plt.savefig(os.path.join(args.out_dir, args.sample + '.png'))


# -- Chunk of coded provided by Luz
# -- overcluster prep. run turbo basic scanpy pipeline
sc.pp.normalize_per_cell(adata,
                         counts_per_cell_after = 1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,
                            min_mean = 0.0125,
                            max_mean = 3,
                            min_disp = 0.5)
# -- Subset to HVGs
adata = adata[:, adata.var['highly_variable']]

sc.pp.scale(adata,
            max_value = 10)

sc.tl.pca(adata,
          svd_solver='arpack')

sc.pp.neighbors(adata)

# -- overclustering proper - do basic clustering first, then cluster each cluster
sc.tl.leiden(adata)
adata.obs['leiden'] = [ str(i) for i in adata.obs['leiden'] ]

for clus in np.unique(adata.obs['leiden']):
    
    adata_sub = adata[adata.obs['leiden'] == clus].copy()
    sc.tl.leiden(adata_sub)
    
    adata_sub.obs['leiden'] = [ clus+','+i for i in adata_sub.obs['leiden'] ]
    adata.obs.loc[adata_sub.obs_names,'leiden'] = adata_sub.obs['leiden']

#compute the cluster scores - the median of Scrublet scores per overclustered cluster
for clus in np.unique(adata.obs['leiden']):
    
    results = np.median(adata.obs.loc[adata.obs['leiden'] == clus, 'scrublet_score'])
    adata.obs.loc[adata.obs['leiden']==clus, 'scrublet_cluster_score'] = results
        
#now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
med = np.median(adata.obs['scrublet_cluster_score'])
mask = adata.obs['scrublet_cluster_score'] > med
mad = np.median(adata.obs['scrublet_cluster_score'][mask]-med)

#let's do a one-sided test. the Bertie write-up does not address this but it makes sense
zscores = (adata.obs['scrublet_cluster_score'].values - med) / (1.4826 * mad)
adata.obs['scrublet_zscore'] = zscores
pvals = 1 - scipy.stats.norm.cdf(zscores)
adata.obs['scrublet_bh_pval'] = kcorr.bh(pvals)
adata.obs['scrublet_bonf_pval'] = kcorr.bonf(pvals)

# -- Extract scrublet annotation
idx = [ bool(re.match('scrublet', i)) for i in adata.obs.columns ]
scrublet_sample = adata.obs.loc[:, idx]

# -- Write results
scrublet_sample.to_csv(os.path.join(args.out_dir, args.sample + 'csv'))