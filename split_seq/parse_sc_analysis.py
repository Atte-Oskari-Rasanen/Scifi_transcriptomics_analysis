import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import os
import scipy.io as sio
import glob
from functools import reduce
import anndata as ad

data_path = '/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/mus/analysis/all-well/DGE_filtered/'

#merge the data from different asyn samples

'''
-------------------------------------------
The class object contains the following properties:
-------------------------------------------
retrieve_data=given the sample path, retrieve the data, do preliminary filtering to reduce array size, return
merge_all=iteratively calls all the dfs using retrieve data, merges them, makes into anndata and returns it

'''
class data_tools():
    def __init__(self, data_path, sample):
        self.data_path=data_path
        self.sample=sample
    def __str__(self):
        return("a class for performing simple processes regarding anndata formatting")

    def get_sample_paths(self):
        samples_list = [dir for dir in sorted(glob.glob(f'{self.data_path}*{self.sample}*')) if os.path.isdir(dir)] 
        return(samples_list)

    def retrieve_data(self, sample_path):
        #sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
        #sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')
        #sc.settings.figdir = '/volume-general/analysis/20201118-novaseq/figures/pbmc/'
        adata1 = sc.read_mtx(sample_path + '/DGE_filtered/DGE.mtx')

        # reading in gene and cell data
        gene_data1 = pd.read_csv(sample_path + '/DGE_filtered/all_genes.csv')
        cell_meta1 = pd.read_csv(sample_path + '/DGE_filtered/cell_metadata.csv')


        # find genes with nan values and filter
        gene_data = gene_data1[gene_data1.gene_name.notnull()]
        notNa = gene_data.index
        notNa = notNa.to_list()

        # remove genes with nan values and assign gene names
        adata = adata1[:,notNa]
        adata.var = gene_data
        adata.var.set_index('gene_name', inplace=True)
        adata.var.index.name = None
        adata.var_names_make_unique()

        # add cell meta data to anndata object
        adata.obs = cell_meta1
        adata.obs.set_index('bc_wells', inplace=True)
        adata.obs.index.name = None
        adata.obs_names_make_unique()
        sc.pp.filter_cells(adata, min_genes=300)
        sc.pp.filter_genes(adata, min_cells=5)

        data1_df=pd.DataFrame.sparse.from_spmatrix(adata.X, index=adata.obs_names, columns=adata.var_names)

        return(data1_df)
    def intersect_genes(self, samples_list):
        all_genes=[]
        for sample_df in samples_list:
            all_genes.append(list(set(sample_df.columns)))
       
        common_genes=set.intersection(*map(set,all_genes))
        # check length
        if len(common_genes) > 0:
            return(common_genes) 
        else:
            return("no common elements")

    def merge_all(self):
        samples_list=self.get_sample_paths()
        all_dfs_sample = list(map(self.retrieve_data, samples_list))
        common_genes = self.intersect_genes(all_dfs_sample)
        df_merged = reduce(lambda df_left,df_right: pd.merge(df_left, df_right,on=common_genes, how='outer'), all_dfs_sample).fillna(0)
        #remove x and y from the non-overlapping ones
        df_merged.drop(df_merged.filter(regex='_y$').columns, axis=1, inplace=True)
        df_merged.drop(df_merged.filter(regex='_x$').columns, axis=1, inplace=True)
        adata = sc.AnnData(df_merged)
        adata.var_names_make_unique()
        return(adata)

a = [1, 2, 3, 4, 5]
b = [5, 6,1, 7, 8,2, 9]
c= [2,5,4]
d=[5,2,4,3]
d=[[1, 2, 3, 4, 5], [5, 6,1, 7, 8,2, 9], [2,5,4],[5,2,4,3]]

t=common_member(a, b)
t
print(common_member(a, b))



data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/analysis_mus/"


dat_format=data_tools(data_path, 'mus')
paths=dat_format.get_sample_paths()
ex_df=dat_format.retrieve_data(paths[2])

all_dfs_sample = list(map(dat_format.retrieve_data, paths))
merged_df=dat_format.merge_all()

def intersect_genes(samples_list):
    all_genes=[]
    for sample_df in samples_list:
        all_genes.append(list(set(sample_df.columns)))
    
    common_genes=set.intersection(*map(set,all_genes))
    # check length
    if len(common_genes) > 0:
        return(common_genes) 
    else:
        return("no common elements")

common_genes = intersect_genes(all_dfs_sample)
len(common_genes)

df_merged

type(all_dfs_sample)
result_1 = pd.concat(all_dfs_sample, join='outer', axis=1).fillna(0)
result_1
# solution 2

df_merged = reduce(lambda df_left,df_right: pd.merge(df_left, df_right,left_index=True, 
right_index=True, how='outer'), all_dfs_sample).fillna(0)



adata=dat_format.merge_all() #warning:AnnData expects .var.index to contain strings, but got values like:[0,1,3]



results_file = '/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/anndata_objects/mus_asyn_sn.h5ad'  # the file that will store the analysis results
adata.write(results_file, compression="gzip")

##############
#Quality control
##############
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata.obs['n_genes_by_counts']
adata.obs['total_counts']

# Scanpy will prepend the string in the save argument with "violin"
# and save it to our figure directory defined in the first step.
'''
the number of genes expressed in the count matrix (so how often a certain gene is expressed in the matrix?)
the total counts per cell
the percentage of counts in mitochondrial genes
'''

sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4) #gene counts is 1800
sc.pl.violin(adata, ['total_counts'], jitter=0.4) #total counts is 4000...
sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4)

#Filtering
# Filter the data
adata = adata[adata.obs.n_genes_by_counts < 800,:]
adata = adata[adata.obs.total_counts < 1500,:]
#adata = adata[adata.obs.pct_counts_mt < 15,:]
adata.shape # Checking number of cells remaining
sc.pl.highest_expr_genes(adata, n_top=20, )

#Visualise
prel_info=[]
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

print('median transcript count per cell: ' + str(adata.obs['tscp_count'].median(0)))
print('median gene count per cell: ' + str(adata.obs['gene_count'].median(0)))
prel_info.append(str(adata.obs['tscp_count'].median(0)))
prel_info.append(str(adata.obs['gene_count'].median(0)))

df = pd.DataFrame(prel_info, columns=['median transcript count per cell','median gene count per cell'])


#Normalise
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
sc.pl.highly_variable_genes(adata, save='') # scanpy generates the filename automatically

# Save raw expression values before variable gene subset
adata.raw = adata

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)


#PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50) # scanpy generates the filename automatically

#UMAP+LEIDEN
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save='_leiden')
