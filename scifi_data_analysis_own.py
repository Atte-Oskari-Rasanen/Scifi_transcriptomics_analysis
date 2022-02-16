#import desc as DESC
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


p = '/media/data/AtteR/scifi-analysis/outputs_starsolo/Scifi5/Scifi_library_2_S2_SciFi5_WP/GeneFull/filtered/'

matrix_filtered = sc.read_mtx(p + "matrix.mtx").T

genes = pd.read_csv(p + 'features.tsv', header=None, sep='\t')
matrix_filtered.var['gene_ids'] = genes[0].values
matrix_filtered.var['gene_symbols'] = genes[1].values
matrix_filtered.var_names = matrix_filtered.var['gene_symbols']
matrix_filtered
# Make sure the gene names are unique
matrix_filtered.var_names_make_unique(join="-")

matrix_filtered_df = matrix_filtered.to_df().head()
matrix_filtered_df.head()

matrix_filtered_df.columns


match_cols = [col for col in matrix_filtered_df.columns if 'AAV_ASYN' in col] #none found
print(match_cols)



cells = pd.read_csv(p + 'barcodes.tsv', header=None, sep='\t')
matrix_filtered.obs['barcode'] = cells[0].values
matrix_filtered.obs_names = cells[0]
# Make sure the cell names are unique
matrix_filtered.obs_names_make_unique(join="-")


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

    matrix_filtered = sc.read_mtx(path + "matrix.mtx").T

    genes = pd.read_csv(path + 'features.tsv', header=None, sep='\t')
    matrix_filtered.var['gene_ids'] = genes[0].values
    matrix_filtered.var['gene_symbols'] = genes[1].values
    matrix_filtered.var_names = matrix_filtered.var['gene_symbols']
    # Make sure the gene names are unique
    matrix_filtered.var_names_make_unique(join="-")

    cells = pd.read_csv(path + 'barcodes.tsv', header=None, sep='\t')
    matrix_filtered.obs['barcode'] = cells[0].values
    matrix_filtered.obs_names = cells[0]
    # Make sure the cell names are unique
    matrix_filtered.obs_names_make_unique(join="-")
    print("imported odt")                          # write a cache file for faster subsequent reading

    #########################3
    matrix_filtered_W = sc.read_mtx(path.replace("oDT", "WP") + "matrix.mtx").T

    genes = pd.read_csv(path.replace("oDT", "WP") + 'features.tsv', header=None, sep='\t')
    matrix_filtered_W.var['gene_ids'] = genes[0].values
    matrix_filtered_W.var['gene_symbols'] = genes[1].values
    matrix_filtered_W.var_names = matrix_filtered_W.var['gene_symbols']
    # Make sure the gene names are unique
    matrix_filtered_W.var_names_make_unique(join="-")

    cells = pd.read_csv(path.replace("oDT", "WP") + 'barcodes.tsv', header=None, sep='\t')
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
matrix_filtered = matrix_filtered.transpose()
matrix_filtered_W

matrix_filtered_W = matrix_filtered_W.transpose()

sc.pl.highest_expr_genes(matrix_filtered_W, n_top=20, )

'''
When creating the matrix afterwards you do the normal gene x cell matrix and add rows in the end corresponding to the viral bcs.
Keep in mind that some of the viral ones can be found in oDT too due to technical reasons. Then you import the viral data separately 
as gene x cell, count the number of them, and import this data to the rows in the original oDT gene x cell  matrix.
'''


#print(matrix_filtered_W['gene_ids'])
print("check if column names of interest exist:")
matrix_filtered_W_df = matrix_filtered_W.to_df().head()
matrix_filtered_df = matrix_filtered.to_df().head()

matrix_filtered_W_df.columns


match_cols = [col for col in matrix_filtered_W_df.columns if 'AAV_ASYN' in col] #none found
print(match_cols)
print("done")

    
#import re                                  # Add the re import declaration to use regex
#reg = re.compile(r'AAV')                    # Compile the regex
#AAVs = ['AAV_CTE', 'AAV_ASYN', 'AAV_BFP', 'AAV_H2BGFP']
# iterating the columns
#index_names = list(filter(lambda x: reg.searc(x), AAVs))
#index_names = list(filter(reg.match, AAVs))
#index_names

matrix_filtered_W_df.index
matrix_filtered_W_df.index.name = 'gene'
matrix_filtered_W_df2=matrix_filtered_W_df
matrix_filtered_W_df2.index.name = 'gene'

matrix_filtered_W_df2
#matrix_filtered_W_df2.reset_index()


matrix_filtered_W_df.reset_index(inplace=True)#print(index_names)
matrix_filtered_W_df['gene']

AAV_data = matrix_filtered_W_df.loc[["AAV_CTE", "AAV_ASYN", "AAV_H2BGFP", "AAV_BFP"]]
AAV_data

#index_names = matrix_filtered_W_df[ (matrix_filtered_W_df['gene'] == 'AAV_CTE') | (matrix_filtered_W_df['gene'] == 'AAV_ASYN') | (matrix_filtered_W_df['gene'] == 'AAV_BFP') | (matrix_filtered_W_df['gene'] == 'AAV_H2BGFP')].index
#index_names

matrix_filtered_df = matrix_filtered.to_df().head()
matrix_filtered_W
matrix_filtered_W_df = matrix_filtered_W_df.loc[index_names]
matrix_filtered_W_df

matrix_filtered_W = matrix_filtered_W

# merge and sum
# we merge the oDT matrix with WP matrix that has been filtered to only contain the AAV genes (rows are genes, cols are bcs)



matrix_filtered_F = matrix_filtered_df.T.merge(AAV_data.T,left_index=True, right_index=True, how='left').T
#matrix_filtered_F = matrix_filtered_df.T.merge(matrix_filtered_W.T,left_index=True, right_index=True, how='left').T
matrix_filtered_F
matrix_filtered_F.index
#setting gene_ids column as index but we dont have it
matrix_filtered_F = matrix_filtered_F.set_index('gene').fillna(0)

matrix_filtered_F.index = matrix_filtered_F.index.str.upper()

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

annMatrix = DF_All.T.iloc[:,1:]
annMatrix
annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : completeTable.T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})
annVar = pd.DataFrame(index=completeTable.iloc[1:,0:1].index, data=completeTable.iloc[1:,0:1].index, columns=['Gene'])
adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata_TX.obs_names_make_unique(join="-")

adata_TX


#get the values in the right other
compT = completeTable.T
compT_in = compT[compT['0_Sample'] == 'SciFi5', :]
compT_in
compT

#filter 
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)



adata = desc.utilities.read_mtx('data/pbmc/matrix.mtx').T

genes = pd.read_csv('data/pbmc/genes.tsv', header=None, sep='\t')
adata.var['gene_ids'] = genes[0].values
adata.var['gene_symbols'] = genes[1].values
adata.var_names = adata.var['gene_symbols']
# Make sure the gene names are unique
adata.var_names_make_unique(join="-")

import scanpy as sc 
for sample_id in samples_list :
    #matrix_filtered_W = matrix_filtered_W.tocsr()

    #matrix_filtered_W = pd.read_pickle(path.replace("oDT", "WP") +"matrix.mtx.gz")
    # Select only AAVs from WPRE cDNA
    #print(type(matrix_filtered_W))
    print(matrix_filtered_W.to_df().head())
    print(matrix_filtered_W)
    
    #matrix_filtered_W.X[:,]
   # print(matrix_filtered_W[:, 'AAV_CTE'].X.tolist())

    #print(matrix_filtered_W['gene_ids'])
    print("check if column names of interest exist:")
    spike_cols = [col for col in matrix_filtered_W.columns if 'AAV' in col]
    print(spike_cols)
    print("done")
    #are the index names only asyn and tagbfp now?
    index_names = matrix_filtered_W[ (matrix_filtered_W['gene_ids'] == 'AAV_CTE') | (matrix_filtered_W['gene_ids'] == 'AAV_ASYN') | (matrix_filtered_W['gene_ids'] == 'AAV_BFP') | (matrix_filtered_W['gene_ids'] == 'AAV_H2BGFP')].index
    #print(index_names)
    matrix_filtered_W = matrix_filtered_W.loc[index_names]
    # merge and sum
    matrix_filtered_F = matrix_filtered.T.merge(matrix_filtered_W.T,left_index=True, right_index=True, how='left').T
    matrix_filtered_F = matrix_filtered_F.set_index('gene_ids').fillna(0)
    matrix_filtered_F.index = matrix_filtered_F.index.str.upper()
    matrix_filtered_F = matrix_filtered_F.groupby(matrix_filtered_F.index).sum()
    #add sampleName
    topRow = pd.DataFrame(columns=matrix_filtered_F.columns, index=['0_Sample']).fillna(sample_Name)
    #topRow.loc['_Sample'] = sample_Name
    matrix_filtered_out = pd.concat([topRow, matrix_filtered_F])
    #Merge into a compete dataframe
    completeTable = completeTable.merge(matrix_filtered_out,left_index=True, right_index=True, how='outer').fillna(0)
    















    #matrix_filtered = pd.read_pickle(path + "matrix.mtx.gz")
    #matrix_filtered = scio.mmread(path + "matrix.mtx")
    #matrix_filtered = sc.read(path + "matrix.mtx")
    #matrix_filtered = sc.read_10x_mtx(
    #    path, # the directory with the `.mtx` file
    #    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    #    cache=True)          

    #matrix_filtered.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`                    # write a cache file for faster subsequent reading
    #matrix_filtered_W = sc.read_10x_mtx(
    #path.replace("oDT", "WP") , # the directory with the `.mtx` file
    #var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    #cache=True)                              # write a cache file for faster subsequent reading
    #matrix_filtered_W.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`                    # write a cache file for faster subsequent reading
