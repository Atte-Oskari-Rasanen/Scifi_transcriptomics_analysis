from turtle import color
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
import anndata as ad
import os
import plotly 
import plotly.graph_objs as go   #numba needs a numpy of 1.21 or less 
import pickle
import glob
import ntpath
import seaborn as sns
import scipy
import scipy.io as scio
import re



basePath = 'SoloNova10x_correctAAV/'
basePath= "/media/data/AtteR/scifi-analysis/outputs_starsolo/Scifi5/"

samples_list = [dir for dir in sorted(glob.glob(basePath +"*")) if os.path.isdir(dir)] 
samples_list= [item for item in samples_list if '_WP' not in item]

samples_list




for sample_id in samples_list:
    sample_id_W = ntpath.basename(sample_id).replace('_oDT','_WP')
    sample_Name = ntpath.basename(sample_id).replace('_oDT','').replace('Scifi_library_','')
    print(sample_id_W)
    print(sample_Name)
    path = sample_id + "/GeneFull/filtered/"
    print(path)
    path_w = path.replace("oDT", "WP")
    print(path_w)
    #matrix_filtered=pd.read_pickle(path + "matrix.mtx")
    #matrix_filtered_W=pd.read_pickle(path.replace("oDT", "WP") + "matrix.mtx")

    matrix_filtered = sc.read_mtx(path + "matrix.mtx.gz").T

    genes = pd.read_csv(path + 'features.tsv.gz', header=None, sep='\t')
    matrix_filtered.var['gene_ids'] = genes[0].values
    matrix_filtered.var['gene_symbols'] = genes[1].values
    matrix_filtered.var_names = matrix_filtered.var['gene_symbols']
    # Make sure the gene names are unique
    matrix_filtered.var_names_make_unique(join="-")

    cells = pd.read_csv(path + 'barcodes.tsv.gz', header=None, sep='\t')
    matrix_filtered.obs['barcode'] = cells[0].values
    matrix_filtered.obs_names = cells[0]
    # Make sure the cell names are unique
    matrix_filtered.obs_names_make_unique(join="-")
    print("imported odt")                          # write a cache file for faster subsequent reading

    #########################3
    matrix_filtered_W = sc.read_mtx(path.replace("oDT", "WP") + "matrix.mtx.gz").T

    genes = pd.read_csv(path.replace("oDT", "WP") + 'features.tsv.gz', header=None, sep='\t')
    matrix_filtered_W.var['gene_ids'] = genes[0].values
    matrix_filtered_W.var['gene_symbols'] = genes[1].values
    matrix_filtered_W.var_names = matrix_filtered_W.var['gene_symbols']
    # Make sure the gene names are unique
    matrix_filtered_W.var_names_make_unique(join="-")

    cells = pd.read_csv(path.replace("oDT", "WP") + 'barcodes.tsv.gz', header=None, sep='\t')
    matrix_filtered_W.obs['barcode'] = cells[0].values
    matrix_filtered_W.obs_names = cells[0]
    # Make sure the cell names are unique
    matrix_filtered_W.obs_names_make_unique(join="-")

    #matrix_filtered = sc.read_10x_mtx(
    #path , # the directory with the `.mtx` file
    #var_names='gene_ids',                # use gene symbols for the variable names (variables-axis index)
    #cache=True)    
    
    #matrix_filtered_W = sc.read_10x_mtx(
    #path_w , # the directory with the `.mtx` file
    #var_names='gene_ids',                # use gene symbols for the variable names (variables-axis index)
    #cache=True)                              # write a cache file for faster subsequent reading
completeTable = pd.DataFrame()

matrix_filtered
matrix_filtered_W

#Atm rows contain cells and columns the genes. We need to turn them around for later indexing based on gene names

matrix_filtered = matrix_filtered.transpose()
matrix_filtered_W = matrix_filtered_W.transpose()

#sc.pl.highest_expr_genes(matrix_filtered_W, n_top=20, )

'''
When creating the matrix afterwards you do the normal gene x cell matrix and add rows in the end corresponding to the viral bcs.
Keep in mind that some of the viral ones can be found in oDT too due to technical reasons. Then you import the viral data separately 
as gene x cell, count the number of them, and import this data to the rows in the original oDT gene x cell  matrix.
'''


#print(matrix_filtered_W['gene_ids'])
print("check if column names of interest exist:")
matrix_filtered_W_df = matrix_filtered_W.to_df()
matrix_filtered_df = matrix_filtered.to_df()

matrix_filtered_W_df.columns


match_cols = [col for col in matrix_filtered_W_df.columns if 'AAV_ASYN' in col] #none found
print(match_cols)
print("done")

AAV_data = matrix_filtered_W_df.loc[["AAV_CTE", "AAV_ASYN", "AAV_H2BGFP", "AAV_BFP"]]
AAV_data

#index_names = matrix_filtered_W_df[ (matrix_filtered_W_df['gene'] == 'AAV_CTE') | (matrix_filtered_W_df['gene'] == 'AAV_ASYN') | (matrix_filtered_W_df['gene'] == 'AAV_BFP') | (matrix_filtered_W_df['gene'] == 'AAV_H2BGFP')].index
#index_names

# merge and sum
# we merge the oDT matrix with WP matrix that has been filtered to only contain the AAV genes (rows are genes, cols are bcs)



matrix_filtered_F = matrix_filtered_df.T.merge(AAV_data.T,left_index=True, right_index=True, how='left').T
#matrix_filtered_F = matrix_filtered_df.T.merge(matrix_filtered_W.T,left_index=True, right_index=True, how='left').T
matrix_filtered_F
matrix_filtered_F.index


#setting gene_ids column as index but we dont have it so skipping it
#matrix_filtered_F = matrix_filtered_F.set_index('gene_symbols').fillna(0)

#matrix_filtered_F.index = matrix_filtered_F.index.str.upper()

#group based on the indeces (genes) and sum
matrix_filtered_F_g = matrix_filtered_F.groupby(matrix_filtered_F.index).sum()
matrix_filtered_F_g


#add sampleName so scifi5, 6 etc.
topRow = pd.DataFrame(columns=matrix_filtered_F_g.columns, index=['0_Sample']).fillna(sample_Name)
topRow 
#topRow.loc['_Sample'] = sample_Name
matrix_filtered_out = pd.concat([topRow, matrix_filtered_F_g])
#Merge into a compete dataframe
completeTable = completeTable.merge(matrix_filtered_out,left_index=True, right_index=True, how='outer').fillna(0)

DF_All = completeTable
DF_All
DF_All.T.index

########################################
########################################


#categorise cells based on the number of a syn copies in them to see marker gene expression
#DF_All = DF_All.drop('0_Sample', axis=0) #remove sample row at least for now as it interferes with generation of anndata object

DF_All.index = map(str.upper, DF_All.index)
DF_All = DF_All.drop('0_SAMPLE', axis=0)


#originally a full table (name was CompleteTable) was created but no need I think
#Fulltable = completeTable.merge(final_df_T_WORKS,left_index=True, right_index=True, how='outer').fillna(0)
#Fulltable

annMatrix = DF_All.T.iloc[:-1,:]
#annMatrix = annMatrix.drop("asyn_copies", axis=1) #.iloc[:,1:]get_level_values('first')
annMatrix
annMatrix.index
#final_df_T_g.T['asyn_copies'].tolist()

#final_df_T_g    # had to add final_df_T_WORKS.T.iloc[:-1,0:1] -1 instead of just all since the last row is rowindeces
annObs = pd.DataFrame(index=DF_All.T.iloc[:-1,0:1].index, data={'CellBarcode':DF_All.T.iloc[:-1,0:1].index})

#annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : .T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})

annObs
annObs.index
annVar = pd.DataFrame(index=DF_All.iloc[:,0:1].index, data=DF_All.iloc[:,0:1].index, columns=['Gene'])
annVar.index
#adata = ad.AnnData(annMatrix, obs=annObs)

#the levels are ogranised differently between annmatrix and annobs
adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata_TX.obs_names_make_unique(join="-")

adata_TX


adata_TX_red = adata_TX.copy()


#knee plot to assess the number of cells we'll be filtering
# Create the "knee plot"
from datetime import datetime 
from datetime import date
run_date = date.today()
str(run_date)

knee = np.sort((np.array(annMatrix.T.sum(axis=1))).flatten())[::-1] #gotta transpose as annmatrix normally has genes as columns but we want cells
fig, ax = plt.subplots(figsize=(10, 7))

ax.loglog(knee, range(len(knee)),linewidth=5, color="g")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
#plt.savefig("/media/data/AtteR/scifi-analysis/Python-scifi-analysis/plots/Kneeplot_"+ str(run_date) + ".png")

plt.show()
plt.close()    # close the figure window

adata_TX_red.obs['total_counts'] = np.sum(adata_TX_red.X, axis=1)  #earlier if was adata_TX
gene_counts = np.sum(adata_TX_red.X, axis=1)   #Take each column (cell) and count the number of transcripts

#How to count the n_genes and n_tot per cell?
adata_TX_red.obs['total_counts']

adata_TX_red
sc.pl.highest_expr_genes(adata_TX_red, n_top=20, )
adata_TX_red


#adata_TX_red.obs['n_genes']  

sc.pp.filter_cells(adata_TX_red, min_genes=200)
sc.pp.filter_genes(adata_TX_red, min_cells=3)
adata_TX_red

adata_TX_red.var['mt'] = adata_TX_red.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_TX_red, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#we would need the number of genes in a cell and the total number of molecules (n_count)
sc.pl.violin(adata_TX_red, ['n_genes', 'total_counts', 'total_counts_mt'],jitter=0.4, multi_panel=True)
#adata_TX_red = adata_TX_red[adata_TX_red.obs['read_counts'] >= 1000].copy()

#adata_TX_red[adata_TX_red.obs['Sample'] == 'SciFi6', :]

# plot histogram of transcript count per cell


plt.savefig("/media/data/AtteR/scifi-analysis/Python-scifi-analysis/plots/Scifi5_violincounts_"+ str(run_date) + ".png")
plt.show()
plt.close()    # close the figure window

#remove cells that contain mt genes 
#sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
#sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# delete cells witl less then 5000 transcripts from visualization
# gene_counts_filt = np.delete(gene_counts, np.where(gene_counts < 5000))
# print(len(gene_counts_filt))
print("Minimum number of transcripts per cell:", np.min(gene_counts), 
      "\n Median number of transcripts per cell:", np.median(gene_counts),
      "\n Maximum number of transcripts per cell:", np.max(gene_counts))
plt.hist(gene_counts, bins=100)
#plt.savefig('SciFiAll_Histogram_2200.pdf')
plt.savefig("/media/data/AtteR/scifi-analysis/Python-scifi-analysis/plots/Scifi5_hist_"+ str(run_date) + ".png")

plt.show()
plt.close()    # close the figure window

# create a backup anndata object 
adata_TX_ref = adata_TX.copy()

#sc.pl.scatter(adata_TX_red, x='total_counts', y='n_genes_by_counts')
adata_TX_red
matplotlib.pyplot.scatter(adata_TX_red.obs['total_counts'], adata_TX_red.obs['n_genes'])
plt.show()
plt.close()
#

#remove outliers 
adata_TX_red = adata_TX_red.obs[adata_TX_red.obs['n_genes'] < 5800, :]
adata_TX_red


#normalise counts
sc.pp.normalize_total(adata_TX_red, target_sum=1e4)
#sc.pp.normalize_per_cell(adata_TX_red, counts_per_cell_after=1e4)
sc.pp.log1p(adata_TX_red)

#find the highly variable genes
sc.pp.highly_variable_genes(adata_TX_red, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata_TX_red)


#filter based on the highly var features 
adata_TX_red_f = adata_TX_red[:, adata_TX_red.var.highly_variable]
sc.pp.regress_out(adata_TX_red_f, ['total_counts', 'pct_counts_mt'])

#Scale each gene to unit variance
sc.pp.scale(adata_TX_red_f, max_value=10)

adata_pca = adata_TX_red_f


#PCA
sc.tl.pca(adata_pca, svd_solver='arpack')
#scatter plot of the pca coordinates
sc.pl.pca(adata_pca, color='n_genes')
adata_pca

#inspect the contribution of the total no of PCAs to the variance in data 
sc.pl.pca_variance_ratio(adata_pca, log=True)
#results_file = "/media/data/AtteR/scifi-analysis/Python-scifi-analysis/scifi5.h5ad"
#adata_pca.write(results_file)


# neighborhood graphs embedding - tried different ones but separation still poor.
sc.pp.neighbors(adata_pca, n_neighbors=20, n_pcs=30)
#
sc.tl.leiden(adata_pca)
# had to add leiden prior to next line for some reason
sc.tl.paga(adata_pca)
sc.pl.paga(adata_pca, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata_pca, init_pos='paga')

#n_componenents is the dims equivalent of seurat
sc.tl.umap(adata_pca, min_dist=0.5, spread=0.1, n_components=2)
sc.tl.umap(adata_pca, n_components=10)

sc.pl.umap(adata_pca, color=['A_SYN'], use_raw=False)

####################


sc.pl.umap(adata_pca)
adata_pca.var["Gene"]
sc.pl.umap(adata_pca, color=['FOXA2', 'CALB1', 'CHST8'])


#As we set the .raw attribute of adata, the previous plots showed the “raw” (normalized, logarithmized, but uncorrected) gene expression. 
# You can also plot the scaled and corrected gene expression by explicitly stating that you don’t want to use .raw.

sc.pl.umap(adata_pca, color=['TH', 'DAT', 'CALB1'], use_raw=False)
sc.pl.umap(adata_pca, color=['leiden'], use_raw=False)

# Leiden clustering
sc.tl.leiden(adata_pca)

sc.pl.umap(adata_pca, color=['leiden'])
sc.pl.umap(adata_pca, color=['AAV_ASYN'])

adata_pca
adata_pca.var['Gene']
#sc.pl.umap(adata_pca, color=['N_asyn'])




#Find marker genes 
sc.tl.rank_genes_groups(adata_pca, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_pca, n_genes=25, sharey=False)


# create a slot with raw counts
adata_pca.raw = adata_pca.copy()
sc.pp.log1p(adata_pca)

adata_pca
sc.pp.scale(adata_pca, max_value=10)
sc.tl.pca(adata_pca, svd_solver='arpack', n_comps=80, use_highly_variable=False)

sc.pl.violin(adata_pca, ['n_genes', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)


# only keep the top 30 PC
adata_pca.obsm['X_pca'] = adata_pca.obsm['X_pca'][:, :30]
adata_pca.varm['PCs'] = adata_pca.varm['PCs'][:, :30]
    
#compute neighbours with top 30
sc.pp.neighbors(adata_pca, n_pcs=20, n_neighbors=40, random_state=1)
sc.tl.tsne(adata_pca, perplexity=20, random_state=1)
sc.tl.umap(adata_pca, min_dist = 0.8, spread = 1.5, n_components=3, random_state=1)

#import os
#os.system("pip3 install leidenalg")


sc.tl.leiden(adata_pca, resolution=1.2, key_added = 'leiden_r0125', random_state=1) # change resolution if desired. the prev val was 0.125
#figsize(5,5)
#sc.pl.tsne(adata_TX_red, color=['leiden_r0125', 'N(A_syn)'], frameon=False, save='allSamples')

sc.pl.umap(adata_pca, color=['leiden_r0125', 'AAV_ASYN', 'TH'], frameon=False)

adata_pca
#sc.pl.tsne(adata_TX_red, color=['Th','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')
sc.pl.tsne(adata_pca, color=['TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')


# if I dont subset the data prior, some of these genes end up being removed probably and thus wont be found in the analysis
sc.pl.umap(adata_pca, color=['TH','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')



#Analysing subgroups 
