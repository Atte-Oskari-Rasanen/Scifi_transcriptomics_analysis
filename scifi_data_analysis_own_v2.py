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

    
#import re                                  # Add the re import declaration to use regex
#reg = re.compile(r'AAV')                    # Compile the regex
#AAVs = ['AAV_CTE', 'AAV_ASYN', 'AAV_BFP', 'AAV_H2BGFP']
# iterating the columns
#index_names = list(filter(lambda x: reg.searc(x), AAVs))
#index_names = list(filter(reg.match, AAVs))
#index_names

#matrix_filtered_W_df.index
#matrix_filtered_df.index
#matrix_filtered_W_df.index.name = 'gene'
#matrix_filtered_W_df2=matrix_filtered_W_df
#matrix_filtered_W_df2.index.name = 'gene'

#matrix_filtered_W_df2
#matrix_filtered_W_df2.reset_index()


#matrix_filtered_W_df.reset_index(inplace=True)#print(index_names)
#matrix_filtered_W_df['gene']
#AAV_data = matrix_filtered_W_df.loc[:, ["AAV_CTE", "AAV_ASYN", "AAV_H2BGFP", "AAV_BFP"]]
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

#SCANPY
########################################
########################################
#make the anndata matrix
annMatrix = DF_All.T.iloc[:,1:]
annMatrix
annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : completeTable.T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})
annObs
annVar = pd.DataFrame(index=completeTable.iloc[1:,0:1].index, data=completeTable.iloc[1:,0:1].index, columns=['Gene'])
annVar
adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata_TX.obs_names_make_unique(join="-")

adata_TX

# remove cells with less than given amount of genes from the main frame if desired
adata_TX_red = adata_TX.copy()
adata_TX_red.obs['read_counts'] = np.sum(adata_TX.X, axis=1)
adata_TX_red = adata_TX_red[adata_TX_red.obs['read_counts'] >= 2200].copy()

adata_TX_red[adata_TX_red.obs['Sample'] == 'SciFi6', :]

# plot histogram of transcript count per cell

gene_counts = np.sum(adata_TX_red.X, axis=1)
gene_counts.shape

# delete cells witl less then 5000 transcripts from visualization
# gene_counts_filt = np.delete(gene_counts, np.where(gene_counts < 5000))
# print(len(gene_counts_filt))
print("Minimum number of transcripts per cell:", np.min(gene_counts), 
      "\n Median number of transcripts per cell:", np.median(gene_counts),
      "\n Maximum number of transcripts per cell:", np.max(gene_counts))
plt.hist(gene_counts, bins=100)
#plt.savefig('SciFiAll_Histogram_2200.pdf')
#plt.show()



########################################
########################################


#categorise cells based on the number of a syn copies in them to see marker gene expression

adata_TX_df= DF_All
adata_TX_df = adata_TX_df.T

DF_All_T = DF_All.T
DF_All_T


cells = DF_All_T.shape[0]
cells
groups_list = []

#Manual approach
#go over each row of the column AAV_ASYN, 
for i in range(cells):
    print(i)
    # printing the third element of the column
    print (DF_All_T['AAV_ASYN'][i])
    groups_list.append(DF_All_T['AAV_ASYN'][i])
    #if i==10:
    #    break

groups_list
len(groups_list)
cells
DF_All_T['asyn_copies'] = groups_list
DF_All_T

#organise so that groups having same numbers are together
DF_All_T_reorg = DF_All_T.sort_values(by = 'asyn_copies')
#go over the row values, compare the value to the prev one, when the value changes, slice the df at this point!

dfs_grouped = {}
i=0
io = 0 #starting index for slicing, the ending one will be defined within the for loop
for i in range(cells-1):
    cell1 = DF_All_T_reorg["asyn_copies"][i]
    cell2 = DF_All_T_reorg["asyn_copies"][i+1]
    
    if cell1!=cell2:
        print(cell1)
        print(cell2)
        asyn_cop = str(int(cell1)) 
        df_name = asyn_cop
        print(df_name)
        dfs_grouped[df_name]=DF_All_T_reorg[io:i]
        io=i  #make 

    #groups_list.append(DF_All_T['AAV_ASYN'][i])

for key in dfs_grouped.keys():
    print(key)

#now we have dfs for each individual group. Should we have instead used the DF_All_T_reorg and just use the copy number groups as equivalents of "samples" used prior?
dfs_grouped["5"]
dfs_grouped["2"]

pd1 = dfs_grouped["1"]
pd2 = dfs_grouped["2"]
pd3 = dfs_grouped["3"]
pd4 = dfs_grouped["4"]
pd5 = dfs_grouped["5"]

final_df = pd.concat([pd1,pd2,pd3,pd4,pd5], axis=0)
final_df_T = final_df.T  #so genes are indeces, cells are cols
final_df_T
fd = final_df.copy()
fd = fd.T


# adding a level to the column, i.e. grouping cells based on the number of asyn copies

groups = final_df['asyn_copies'].tolist()
len(groups)

#generate a multi index for the df and assign it
final_df_T.columns = pd.MultiIndex.from_arrays([groups, final_df_T.columns])
final_df_T.columns.levels[0]
final_df_T.head()

#final_df2=final_df_T.T
#final_df2.head()

#Other tattemps, didnt work though
groups_df = pd.MultiIndex.from_product([col_names], names =[groups])
groups_df
groups_df = pd.MultiIndex.from_product(groups, col_names)

groups = pd.DataFrame(columns=final_df.columns, index=['asyn_copies'])
groups 
#topRow.loc['_Sample'] = sample_Name
matrix_filtered_out = pd.concat([topRow, matrix_filtered_F_g])

groups = final_df['asyn_copies'].tolist()
groups

topRow = pd.DataFrame(columns=final_df.columns, index=['asyn_copies']).fillna(sample_Name)


final_df = pd.concat(dfs_grouped["0"], dfs_grouped["1", dfs_grouped["2"], dfs_grouped["3"], dfs_grouped["4"], dfs_grouped["5"]], axis=0)

final_df = pd.concat([dfs_grouped["0"], dfs_grouped["1", dfs_grouped["2"], dfs_grouped["3"], dfs_grouped["4"], dfs_grouped["5"]]], axis=0)

final_df = pd.concat([dfs_grouped["0"], dfs_grouped["1", dfs_grouped["2"], dfs_grouped["3"], dfs_grouped["4"], dfs_grouped["5"]]], ignore_index=True)



topRow = pd.DataFrame(columns=matrix_filtered_F_g.columns, index=['0_Sample']).fillna(sample_Name)


###############











# create a backup anndata object 
adata_TX_ref = adata_TX_red.copy()

# create a slot with raw counts
adata_TX_red.raw = adata_TX_red.copy()
sc.pp.log1p(adata_TX_red)

sc.pp.scale(adata_TX_red, max_value=10)
sc.tl.pca(adata_TX_red, svd_solver='arpack', n_comps=80, use_highly_variable=False)

# only keep the top 30 PC
adata_TX_red.obsm['X_pca'] = adata_TX_red.obsm['X_pca'][:, :30]
adata_TX_red.varm['PCs'] = adata_TX_red.varm['PCs'][:, :30]
    
#compute neighbours with top 30
sc.pp.neighbors(adata_TX_red, n_pcs=30, n_neighbors=40, random_state=1)
sc.tl.tsne(adata_TX_red, perplexity=20, random_state=1)
sc.tl.umap(adata_TX_red, min_dist = 0.8, spread = 1.5, n_components=3, random_state=1)

import os
os.system("pip3 install leidenalg")


sc.tl.leiden(adata_TX_red, resolution=0.125, key_added = 'leiden_r0125', random_state=1) # change resolution if desired, this value is low
figsize(5,5)
sc.pl.tsne(adata_TX_red, color=['leiden_r0125', 'Sample'], frameon=False, save='allSamples')

sc.pl.umap(adata_TX_red, color=['leiden_r0125', 'Sample'], frameon=False, save='allSamples')

sc.pl.tsne(adata_TX_red, color=['TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')

sc.pl.umap(adata_TX_red, color=['Th','Slc6a3','Kcnj6','Calb1', 'Nr4a2', 'Drd2','Gfap','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')



#Analysing subgroups 

annMatrix = DF_All.T.iloc[:,1:]
annMatrix
annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : completeTable.T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})
annVar = pd.DataFrame(index=completeTable.iloc[1:,0:1].index, data=completeTable.iloc[1:,0:1].index, columns=['Gene'])
adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata_TX.obs_names_make_unique(join="-")

adata_TX
