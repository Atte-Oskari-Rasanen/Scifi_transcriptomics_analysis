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
DF_All.T.index

########################################
########################################


#categorise cells based on the number of a syn copies in them to see marker gene expression
#DF_All = DF_All.drop('0_Sample', axis=0) #remove sample row at least for now as it interferes with generation of anndata object
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

final_df_T.columns

#the issue I had was due to the multiindexing and thus generating columns with multiple layers. earlier the sample column was a row of values instead. gotta do the same.
#generate a multi index for the df and assign it


# A functional approach

groups = final_df_T.iloc[-1,:] #last row, i.e. asyn copies and all cols

#groups2 = groups.drop("0_Sample", axis=1)
#Series is a type of list in pandas which can take integer values, string values, double values and more. groups data type was a list initially so had to transform into df!!
#But in Pandas Series we return an object in the form of list, having index starting from 0 to n, Where n is the length of values in series.
# series can take a single list, i.e. it is a list but a pd version of it. df can take multiple series.
groups_df = pd.DataFrame(groups).T
groups_df.columns

groups_l = groups_df.values.tolist()

final_df_T.iloc[:,1] #
final_df_T

# the original format used in the example was:
#
#                     0_Sample A1bg A1cf  A2m A2ml1 A3galt2 A4galt A4gnt AA926063  ... hist1h2ail2 mageb1l1 mrpl11 mrpl24 mrpl9 mt-Rnr1 mt-Rnr2 rnf141 tGap1
#AAAGATGAGACGAAAG  2_S2_SciFi5  0.0  0.0  0.0   0.0     0.0    0.0   0.0      0.0  ...         0.0      0.0    0.0    1.0   1.0     0.0     1.0    0.0   0.0
#.... 
#
# replace 0_sample with the asyn copy number

final_df_T.index


#def shit_row_to_top(df, index_to_shift):
#    idx = [i for i in df.index if i!=index_to_shift]
#    return df.loc[idx]

#groups_df
#groups_df1 = groups_df.drop('asyn_copies', axis=0) #remove the asyn_copies column, earlier axis was 1
#groups_df1

#final_df_T_WORKS = pd.concat([final_df_T, groups_df], axis=0)
#final_df_T_WORKS

#extra
# test = final_df_T_WORKS.drop('0_Sample', axis=0)
# test
# #move asyn copies as the first row (right after the title row)
# test["row_indeces"] = range(1,len(test)+1)
# test
# test.iloc[34326,-1] = 0
# test.sort_values("row_indeces").drop('row_indeces', axis=1)
# test
#


######
final_df_T_WORKS = final_df_T.drop('0_Sample', axis=0)
#move asyn copies as the first row (right after the title row)
final_df_T_WORKS["row_indeces"] = range(1,len(final_df_T_WORKS)+1)
final_df_T_WORKS
final_df_T_WORKS.iloc[34326,-1] = 0
final_df_T_WORKS = final_df_T_WORKS.sort_values("row_indeces").drop('row_indeces', axis=1)
final_df_T_WORKS






DA_markers = ['ASYN_COPIES', 'A1BG', 'ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4','GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B','HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8','IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2', 'OTX1', 'SYN2', 'CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83','NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A']
final_df_T_WORKS.index = map(str.upper, final_df_T_WORKS.index)



genes = final_df_T_WORKS.shape[0]
final_df_T_WORKS.index = map(str.upper, final_df_T_WORKS.index)

def subset_df(df, markers):
    df.index = map(str.upper, df.index)
    genes = df.shape[0]

    #remove any duplicates entered into the markers list by accident
    markers = list( dict.fromkeys(markers) )
    markers = list( dict.fromkeys(markers) ) 

    genes = final_df_T_WORKS.shape[0]
    final_df_T_WORKS.index = map(str.upper, final_df_T_WORKS.index)

    subset_dict = dict()
    i = 0
    for gene in range(genes):
        #print(gene)
        if df.index[gene] in markers:
            subset_dict[df.index[gene]] = df.iloc[gene,:]
            i += 1
    print("Found in total " + str(i) + " matches!")
    df_DA = pd.DataFrame.from_dict(subset_dict).T

    return df_DA

df_DA = subset_df(final_df_T_WORKS, DA_markers)
df_DA


'''
Find the genes that are in the list, save these into dic where the key is the gene and value is the cells with the values 
   c c  c c c c 
g
g
g
g
g
'''

#final_df_T

#first need to go over the df, find rownames that do match with da markers, save the info of the matching markers into a list and then
#subset
final_df_T_WORKS

DA_markers = ['A1BG', 'ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4','GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B','HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8','IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2', 'OTX1', 'SYN2', 'CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83','NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A']

final_df_T_WORKS.index = map(str.upper, final_df_T_WORKS.index)
final_df_T_WORKS.index
final_df_T_WORKS_DA = final_df_T_WORKS.loc[DA_markers]
final_df_T_WORKS_DA = final_df_T_WORKS[final_df_T_WORKS[:].isin(DA_markers)]
final_df_T_WORKS_DA

final_df_T_WORKS[:]
final_df_T_WORKS
fda[DA_markers].any()
keyrows = pool[pool.Name.str.contains(key)]


#originally a full table (name was CompleteTable) was created but no need I think
#Fulltable = completeTable.merge(final_df_T_WORKS,left_index=True, right_index=True, how='outer').fillna(0)
#Fulltable

annMatrix = final_df_T_WORKS.T.iloc[:-1,:]
#annMatrix = annMatrix.drop("asyn_copies", axis=1) #.iloc[:,1:]get_level_values('first')
annMatrix
annMatrix.index
#final_df_T_g.T['asyn_copies'].tolist()

#final_df_T_g    # had to add final_df_T_WORKS.T.iloc[:-1,0:1] -1 instead of just all since the last row is rowindeces
annObs = pd.DataFrame(index=final_df_T_WORKS.T.iloc[:-1,0:1].index, data={'CellBarcode':final_df_T_WORKS.T.iloc[:-1,0:1].index, 'N_asyn':final_df_T_WORKS.loc['asyn_copies']})

#annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : .T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})

annObs
annObs.index
annVar = pd.DataFrame(index=final_df_T_WORKS.iloc[:,0:1].index, data=final_df_T_WORKS.iloc[:,0:1].index, columns=['Gene'])
annVar.index
#adata = ad.AnnData(annMatrix, obs=annObs)

#the levels are ogranised differently between annmatrix and annobs
adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata_TX.obs_names_make_unique(join="-")

adata_TX

#subset based on DA_markers 
adata_TX_DA = adata_TX[adata_TX.var['Gene'].isin(['2', '3'])]

adata_TX_red = adata_TX.copy()

#in case you want to look at specific groups:
#adata_TX_red[adata_TX_red.obs['N_asyn'] == '1', :]

DA_markers = ['ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4','GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B','HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8','IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2', 'OTX1', 'SYN2', 'CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83','NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A']

keyrows = pool[pool.Name.str.contains(key)]
# selecting rows based on condition
rslt_df = adata_TX_red[adata_TX_red['Stream'].isin(options)]
#adata_TX_red.obs[:]

adata_TX_red.var[adata_TX_red.var.Gene.isin(DA_markers)]

adata_TX_red.var['Gene']
adata_TX_red.var.Gene.isin(DA_markers)


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
plt.close(fig)    # close the figure window


adata_TX_red.obs['read_counts'] = np.sum(adata_TX.X, axis=1)
gene_counts = np.sum(adata_TX_red.X, axis=1)   #Take each column (cell) and count the number of transcripts


adata_TX_red.obs['N_asyn']
sc.pl.highest_expr_genes(adata_TX_red, n_top=20, )

adata_TX_red
adata_TX_red.obs['n_genes']  

sc.pp.filter_cells(adata_TX_red, min_genes=1000)
sc.pp.filter_genes(adata_TX_red, min_cells=3)
adata_TX_red

adata_TX_red.var['mt'] = adata_TX_red.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_TX_red, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#we would need the number of genes in a cell and the total number of molecules (n_count)
sc.pl.violin(adata_TX_red, ['n_genes'],jitter=0.4, multi_panel=True)
#adata_TX_red = adata_TX_red[adata_TX_red.obs['read_counts'] >= 1000].copy()

#adata_TX_red[adata_TX_red.obs['Sample'] == 'SciFi6', :]

# plot histogram of transcript count per cell

#useless since we already counted these into the obs var as n_genes...
gene_counts = np.sum(adata_TX_red.X, axis=1)   #Take each column (cell) and count the number of transcripts
gene_counts.shape

adata_TX_red2 = adata_TX_red.copy()
adata_TX_red.obs['n_count'] = np.sum(adata_TX_red.X, axis=1)


#from scipy.sparse import csr_matrix
#adata_TX_red.X = csr_matrix(adata_TX_red.X)
#adata_TX_red.write(filename="/media/data/AtteR/scifi-analysis/R-scifi-analysis/objects/Scifi5_anndata.h5ad")

adata_TX_red

sc.pl.violin(adata_TX_red, ['n_genes', 'n_count'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata_TX_red, ['n_genes', 'n_count', 'pct_counts_mt'],jitter=0.4, multi_panel=True)
plt.savefig("/media/data/AtteR/scifi-analysis/Python-scifi-analysis/plots/Scifi5_violincounts_"+ str(run_date) + ".png")
plt.show()
plt.close(fig)    # close the figure window

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
matplotlib.pyplot.scatter(adata_TX_red.obs['n_count'], adata_TX_red.obs['n_genes'])
plt.show()
plt.close()
#
adata_TX_red = adata_TX_red.obs[adata_TX_red.obs['n_genes'] < 5800, :]
adata_TX_red


#normalise counts
sc.pp.normalize_per_cell(adata_TX, counts_per_cell_after=1e4)
sc.pp.log1p(adata_TX)


# create a slot with raw counts
adata_TX_red.raw = adata_TX_red.copy()
sc.pp.log1p(adata_TX_red)

sc.pp.scale(adata_TX_red, max_value=10)
sc.tl.pca(adata_TX_red, svd_solver='arpack', n_comps=80, use_highly_variable=False)

sc.pl.violin(adata_TX_red, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)


# only keep the top 30 PC
adata_TX_red.obsm['X_pca'] = adata_TX_red.obsm['X_pca'][:, :30]
adata_TX_red.varm['PCs'] = adata_TX_red.varm['PCs'][:, :30]
    
#compute neighbours with top 30
sc.pp.neighbors(adata_TX_red, n_pcs=30, n_neighbors=40, random_state=1)
sc.tl.tsne(adata_TX_red, perplexity=20, random_state=1)
sc.tl.umap(adata_TX_red, min_dist = 0.8, spread = 1.5, n_components=3, random_state=1)

#import os
#os.system("pip3 install leidenalg")


sc.tl.leiden(adata_TX_red, resolution=0.75, key_added = 'leiden_r0125', random_state=1) # change resolution if desired. the prev val was 0.125
figsize(5,5)
sc.pl.tsne(adata_TX_red, color=['leiden_r0125', 'N(A_syn)'], frameon=False, save='allSamples')
adata_TX_red
sc.pl.umap(adata_TX_red, color=['leiden_r0125', 'N(A_syn)'], frameon=False, save='allSamples')

#sc.pl.tsne(adata_TX_red, color=['Th','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')

sc.pl.umap(adata_TX_red, color=['Th','Slc6a3','Kcnj6','Calb1', 'Nr4a2', 'Drd2','Gfap','AAV_ASYN','AAV_H2BGFP'], color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False, save='_individualGenes')



#Analysing subgroups 
