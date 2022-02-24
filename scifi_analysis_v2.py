

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

DF_All_T=DF_All.T

#earlier DF_All was DF_All_T
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
DA_markers = ['TH', 'DAT', 'A1BG', 'ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4','GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B','HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8','IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2', 'OTX1', 'SYN2', 'CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83','NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A']


pd0 = dfs_grouped["0"]
pd1 = dfs_grouped["1"]
pd2 = dfs_grouped["2"]
pd3 = dfs_grouped["3"]
pd4 = dfs_grouped["4"]
pd5 = dfs_grouped["5"]
a = len(dfs_grouped["0"]) + len(dfs_grouped["1"]) + len(dfs_grouped["2"]) + len(dfs_grouped["3"]) + len(dfs_grouped["4"]) + len(dfs_grouped["5"])
a

final_df = pd.concat([pd0,pd1,pd2,pd3,pd4,pd5], axis=0)
final_df_T = final_df.T  #so genes are indeces, cells are cols
final_df_T.index
final_df_T

fd = final_df.copy()
fd = fd.T


######
final_df_T_WORKS = final_df_T.drop('0_Sample', axis=0)
final_df_T_WORKS = final_df_T
#move asyn copies as the first row (right after the title row)
final_df_T["row_indeces"] = range(1,len(final_df_T_WORKS)+1)
final_df_T_WORKS
final_df_T_WORKS.iloc[34326,-1] = 0
final_df_T_WORKS = final_df_T_WORKS.sort_values("row_indeces").drop('row_indeces', axis=1)
final_df_T_WORKS


#use this when grouped by asyn but no filtering with regards to DA
final_df_T_WORKS
final_df_T_WORKS.index = map(str.upper, final_df_T_WORKS.index)

#use this when raw
DF_All 

annMatrix = final_df_T_WORKS.T.iloc[:-1,:]
#annMatrix = annMatrix.drop("asyn_copies", axis=1) #.iloc[:,1:]get_level_values('first')
annMatrix
annMatrix.index
#final_df_T_g.T['asyn_copies'].tolist()

#final_df_T_g    # had to add final_df_T_WORKS.T.iloc[:-1,0:1] -1 instead of just all since the last row is rowindeces
annObs = pd.DataFrame(index=final_df_T_WORKS.T.iloc[:-1,0:1].index, data={'CellBarcode':final_df_T_WORKS.T.iloc[:-1,0:1].index, 'N_asyn':final_df_T_WORKS.loc['ASYN_COPIES']})

#annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : .T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})

annObs
annObs.index
annVar = pd.DataFrame(index=final_df_T_WORKS.iloc[:,0:1].index, data=final_df_T_WORKS.iloc[:,0:1].index, columns=['Gene'])
annVar.index
#adata = ad.AnnData(annMatrix, obs=annObs)

#the levels are ogranised differently between annmatrix and annobs
adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata_TX.obs_names_make_unique(join="-")
#adata_TX.write(filename="/media/data/AtteR/scifi-analysis/Python-scifi-analysis/plots_no_groups/All.h5ad", compression=None, compression_opts=None, force_dense=None, as_dense=())

adata_TX

adata_TX.obs['CellBarcode']


adata_TX_red = adata_TX.copy()

# mitochondrial genes
adata_TX_red.var['mt'] = adata_TX_red.var_names.str.startswith('MT-') 
# ribosomal genes
adata_TX_red.var['ribo'] = adata_TX_red.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
adata_TX_red.var['hb'] = adata_TX_red.var_names.str.contains(("^HB[^(P)]"))

#sc.pp.calculate_qc_metrics(adata_TX_red, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)

# remove cells with less than given amount of genes from the main frame if desired
adata_TX_red = adata_TX.copy()
adata_TX_red.obs['read_counts'] = np.sum(adata_TX_red.X, axis=1)
adata_TX_red = adata_TX_red[adata_TX_red.obs['read_counts'] >= 2200].copy()
adata_TX_red

gene_counts = np.sum(adata_TX_red.X, axis=1)
gene_counts.shape

# delete cells witl less then 5000 transcripts from visualization
# gene_counts_filt = np.delete(gene_counts, np.where(gene_counts < 5000))
# print(len(gene_counts_filt))
print("Minimum number of transcripts per cell:", np.min(gene_counts), 
      "\n Median number of transcripts per cell:", np.median(gene_counts),
      "\n Maximum number of transcripts per cell:", np.max(gene_counts))
plt.hist(gene_counts, bins=100)
plt.savefig('SciFiDA_Histogram_2200.pdf')
plt.show()

# create a backup anndata object 
adata_TX_ref = adata_TX_red.copy()


# create a slot with raw counts
adata_TX_red.raw = adata_TX_red.copy()
sc.pp.log1p(adata_TX_red)

# regress out unwanted variables
sc.pp.regress_out(adata_TX_red, ['total_counts'])

sc.pp.scale(adata_TX_red, max_value=10)
sc.tl.pca(adata_TX_red, svd_solver='arpack', n_comps=80, use_highly_variable=False)

# only keep the top 30 PC
adata_TX_red.obsm['X_pca'] = adata_TX_red.obsm['X_pca'][:, :30]
adata_TX_red.varm['PCs'] = adata_TX_red.varm['PCs'][:, :30]
adata_TX_red
#compute neighbours with top 30
sc.pp.neighbors(adata_TX_red, n_pcs=20, n_neighbors=40, random_state=1)
sc.tl.umap(adata_TX_red, min_dist = 0.8, spread = 1.5, n_components=3, random_state=1)
adata_TX_red

adata_TX_red.obs['N_asyn'] = adata_TX_red.obs['N_asyn'].astype('category')

sc.pl.violin(adata_TX_red, ["read_counts"], jitter=0.4, groupby = 'N_asyn', rotation= 45)
sc.tl.leiden(adata_TX_red, resolution=1.2, key_added = 'leiden_r12', random_state=1) # change resolution if desired, this value is low
figsize(5,5)
sc.pl.umap(adata_TX_red, color=['leiden_r12', "N_asyn"], frameon=False, save='allSamples')

sc.pl.umap(adata_TX_red, color='N_asyn')

sc.pl.umap(adata_TX_red, color=['TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')
sc.pl.umap(adata_TX_red, color=DA_markers10, color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')
sc.pl.umap(adata_TX_red, color=DA_markers2, color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')

#run with 10 components, save to a new object so that the umap with 2D is not overwritten.
umap10 = sc.tl.umap(adata_TX_red, n_components=10, copy=True)
fig, axs = plt.subplots(1, 3, figsize=(10,4),constrained_layout=True)

sc.pl.umap(adata_TX_red, color='N_asyn',  title="UMAP")
sc.pl.umap(umap10, color='N_asyn', title="UMAP10", show=False, ax=axs[1], components=[1,2])
sc.pl.umap(umap10, color='N_asyn', title="UMAP10", show=False, ax=axs[2], components=[3,4])
# we can also plot the umap with neighbor edges
sc.pl.umap(adata_TX_red, color='N_asyn', title="UMAP", edges=True)

#lets plot the scaled and normalised data
# compute variable genes  0.0125
sc.pp.highly_variable_genes(adata_TX_red, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes: %d"%sum(adata_TX_red.var.highly_variable))


length = len(DA_markers)

middle_index = length//2

DA_markers_1st = DA_markers[:middle_index]

DA_markers_2nd = DA_markers[middle_index:]

sc.pl.pca_loadings(adata_TX_red, components=[1,2,3,4,5,6,7,8])
sc.pl.pca_loadings(adata_TX_red, components=[3,4])

#plot variable genes
sc.pl.highly_variable_genes(adata_TX_red)

# subset for variable genes in the dataset
adata_TX_red = adata_TX_red[:, adata_TX_red.var['highly_variable']]
DA_markers
var_genes = adata_TX_red.var.highly_variable
var_genes.index[var_genes]
varg = [x for x in DA_markers_1st if x in var_genes.index[var_genes]]
sc.pl.umap(adata_TX_red, color=varg, use_raw=False)



########
# calculate cluster enriched genes for leiden clustering
#Different clustering methods
sc.pl.dendrogram(adata_TX_red, groupby = "leiden_r12")
sc.pl.dotplot(adata_TX_red, DA_markers, groupby='leiden_r12', dendrogram=True)


#hierarch clustering
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

# extract pca coordinates
X_pca = adata_TX_red.obsm['X_pca'] 
X_pca

# kmeans with k=5
kmeans = KMeans(n_clusters=5, random_state=0).fit(X_pca) 
adata_TX_red.obs['kmeans5'] = kmeans.labels_.astype(str)

# kmeans with k=10
kmeans = KMeans(n_clusters=10, random_state=0).fit(X_pca) 
adata_TX_red.obs['kmeans10'] = kmeans.labels_.astype(str)

# kmeans with k=15
kmeans = KMeans(n_clusters=15, random_state=0).fit(X_pca) 
adata_TX_red.obs['kmeans15'] = kmeans.labels_.astype(str)

sc.pl.umap(adata_TX_red, color=['kmeans5', 'kmeans10', 'kmeans15'])
