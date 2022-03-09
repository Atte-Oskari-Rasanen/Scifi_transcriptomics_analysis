

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
basePath= "/media/data/AtteR/scifi-analysis/scifi6/2nd_try/starsolo_outputs/"

def import_sc(sample_path):
    sample_id_W = ntpath.basename(sample_path).replace('_oDT','_WP')
    sample_Name = ntpath.basename(sample_path).replace('_oDT','').replace('Scifi_library_','')
    print(sample_Name)
    path = sample_path + "/Solo.out/GeneFull/filtered/"
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

    #########################
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

    matrices = [matrix_filtered, matrix_filtered_W]
    return(matrices)
    #matrix_filtered = sc.read_10x_mtx(

#RT types
#H+: F8, G8, odT_5, WP_5
#H-: H8, A9, odT_6, odT_6

#So we need to compare pools between each other and stuff within them

#start by taking a single pool and group cells based on wheter they are H+ or H-


samples_list = [dir for dir in sorted(glob.glob(basePath +"*")) if os.path.isdir(dir)] 
samples_list= [item for item in samples_list if '_WP' not in item]


P1_Hpos = ["F8", "G8"]
P1_Hneg = ["H8", "A9"]
P2_Hpos = ["oDT_5"]
P2_Hneg = ["oDT_6"]

P_RT_groups = [P1_Hpos, P1_Hneg, P2_Hpos, P2_Hneg]
P_RT_groups


# go over the groups, when the match is found, import the data based on import_sc function with oDT and WP pairs. 
# the list element name (i.e. the group) as the key and imported sc objects (oDT, WP) as value

grouped_sc = {}
for group in P_RT_groups: 
    for d in samples_list:
        f = ntpath.basename(d)
        #print(f)
        check_match = [ele for ele in group if(ele in f)]
        #check_match = [f for f in f if(f in P1_Hpos)]

        #match = any(P2_Hpos in d for P2_Hpos in P2_Hpos)
        if bool(check_match)==True:
            group_name=check_match[0]
            print("True: " + f)
            print(d)

            grouped_sc[group_name]=import_sc(d)


#we have the dicts consisting of RT - df pairs.

#function that takes in the dataframes of oDT and WP and merges them based on our needs
def full_df(oDT_df, WP_df):
    #get the aav vector data   --- do we need cte and bfp tho
    AAV_data = WP_df.loc[["AAV_CTE", "AAV_ASYN", "AAV_H2BGFP", "AAV_BFP"]]
    AAV_data

    #left_index: If True, use the index (row labels) from the left DataFrame or Series as its join key(s)
    #we merge based on gene names meaning that need to transpose the dfs 
    df_full = pd.merge(oDT_df.T, AAV_data.T, left_index=True, right_index=True, how='left').fillna(0)
    #AAV_SYN_x and AAV_SYN_y are generated since there is a clash of columns which were not involved in 
    # the merge operation initially! Need to sum these columns, thus generating a new column which I name as 
    # the original and then remove the x and y suffix containing ones

    ASYN_sum = df_full.loc[:,"AAV_ASYN_x"] + df_full.loc[:,"AAV_ASYN_y"] 
    CTE_sum = df_full.loc[:,"AAV_CTE_x"] + df_full.loc[:,"AAV_CTE_y"]
    H2BGFP_sum = df_full.loc[:,"AAV_H2BGFP_x"] + df_full.loc[:,"AAV_H2BGFP_y"]
    BFP_sum = df_full.loc[:,"AAV_BFP_x"] + df_full.loc[:,"AAV_BFP_x"]

    df_full["AAV_ASYN"] = ASYN_sum
    df_full["AAV_CTE"]= CTE_sum
    df_full["AAV_H2BGFP"] = H2BGFP_sum
    df_full["AAV_BFP"]= BFP_sum
    df_full
    del df_full["AAV_ASYN_x"]
    del df_full["AAV_ASYN_y"]
    del df_full["AAV_CTE_x"]
    del df_full["AAV_CTE_y"]
    del df_full["AAV_H2BGFP_x"]
    del df_full["AAV_H2BGFP_y"]
    del df_full["AAV_BFP_x"]
    del df_full["AAV_BFP_y"]

    df_full = df_full.T
    df_full = df_full.drop_duplicates()
    return(df_full)


#Pool1, RT H+
matrix_filtered_p1_1 = grouped_sc['F8'][0].transpose() #0 is odT, 1 is WP seq
matrix_filtered_W_p1_1 = grouped_sc['F8'][1].transpose()

matrix_filtered_W_df_p1_1_W  = matrix_filtered_W_p1_1.to_df()
matrix_filtered_W_df_p1_1_W
matrix_filtered_df_p1_1  = matrix_filtered_p1_1.to_df()



#index_names = matrix_filtered_W_df_p1_1_W[ (matrix_filtered_W_df_p1_1_W.iloc[0,:]== 'AAV_CTE') | (matrix_filtered_W_df_p1_1_W.iloc[0,:]== 'AAV_ASYN') | (matrix_filtered_W_df_p1_1_W.iloc[0,:] == 'AAV_BFP') | (matrix_filtered_W_df_p1_1_W.iloc[0,:] == 'AAV_H2BGFP')].index

#index_names = matrix_filtered_W_p1_1[ (matrix_filtered_W_df_p1_1['gene_symbols']== 'AAV_CTE') | (matrix_filtered_W_df_p1_1['gene_ids'] == 'AAV_ASYN') | (matrix_filtered_W_df_p1_1['gene_ids'] == 'AAV_BFP') | (matrix_filtered_W_df_p1_1['gene_ids'] == 'AAV_H2BGFP')].index
#print(index_names)
#matrix_filtered_W = matrix_filtered_W_df_p1_1_W.T.loc[index_names]
matrix_filtered_W_df_p1_1_W

matrix_filtered_p1_2 = grouped_sc['G8'][0].transpose()
matrix_filtered_W_p1_2  = grouped_sc['G8'][1].transpose()
matrix_filtered_W_df_p1_2  = matrix_filtered_W_p1_2.to_df()
matrix_filtered_df_p1_2  = matrix_filtered_p1_2.to_df()

#Combine the dataframes 
Pool1_RT_Hpos_df = pd.concat([matrix_filtered_df_p1_1, matrix_filtered_df_p1_2], axis=1)
Pool1_RT_Hpos_df_WP = pd.concat([matrix_filtered_W_df_p1_1_W, matrix_filtered_W_df_p1_2], axis=1)

Pool1_final_RT_Hpos_full = full_df(Pool1_RT_Hpos_df,Pool1_RT_Hpos_df_WP)

#############################################################################
#############################################################################
#############################################################################

#Pool1, RT H-
matrix_filtered_p1_1 = grouped_sc['H8'][0].transpose()
matrix_filtered_W_p1_1 = grouped_sc['H8'][1].transpose()

matrix_filtered_W_df_p1_1  = matrix_filtered_W_p1_1.to_df()
matrix_filtered_df_p1_1  = matrix_filtered_p1_1.to_df()


matrix_filtered_p1_2 = grouped_sc['A9'][0].transpose()
matrix_filtered_W_p1_2  = grouped_sc['A9'][1].transpose()
matrix_filtered_W_df_p1_2  = matrix_filtered_W_p1_2.to_df()
matrix_filtered_df_p1_2  = matrix_filtered_p1_2.to_df()

Pool1_RT_Hneg_df = pd.concat([matrix_filtered_df_p1_1, matrix_filtered_df_p1_2], axis=1)
Pool1_RT_Hneg_df_WP = pd.concat([matrix_filtered_W_df_p1_1, matrix_filtered_W_df_p1_2], axis=1)

Pool1_final_RT_Hneg_full = full_df(Pool1_RT_Hneg_df,Pool1_RT_Hneg_df_WP)



##################
#Pool2, RT H+
Pool2_RT_Hpos = grouped_sc['oDT_5'][0].transpose()
Pool2_RT_Hpos_W = grouped_sc['oDT_5'][1].transpose()

Pool2_RT_Hpos_W_df = Pool2_RT_Hpos_W.to_df()
Pool2_RT_Hpos_df  = Pool2_RT_Hpos.to_df()

Pool2_RT_Hpos_full = full_df(Pool2_RT_Hpos_df,Pool2_RT_Hpos_W_df)



#Pool2, RT H-
Pool2_RT_Hneg = grouped_sc['oDT_6'][0].transpose()
Pool2_RT_Hneg_W = grouped_sc['oDT_6'][1].transpose()

Pool2_RT_Hneg_W_df = Pool2_RT_Hneg_W.to_df()
Pool2_RT_Hneg_df  = Pool2_RT_Hneg.to_df()

Pool2_final_RT_Hneg_full = full_df(Pool2_RT_Hneg_df,Pool2_RT_Hneg_W_df)



'''
##########################################################
#used drop duplicates prior since found plenty of duplicates:
duplicated_columns_list = []
list_of_all_columns = list(Pool1_final_RT_Hpos_full.columns)
for column in list_of_all_columns:
    if list_of_all_columns.count(column) > 1 and not column in duplicated_columns_list:
        duplicated_columns_list.append(column)
duplicated_columns_list
Pool1_final_RT_Hpos_full['TGTACCCTGGGTGCAA']  #duplicates
##########################################################
'''
#Generate anndata matrices
def df_to_anndata(df):
    annMatrix = df.T.iloc[:-1,:]
    #annMatrix = annMatrix.drop("asyn_copies", axis=1) #.iloc[:,1:]get_level_values('first')
    annMatrix
    annMatrix.index
    #final_df_T_g.T['asyn_copies'].tolist()

    #final_df_T_g    # had to add final_df_T_WORKS.T.iloc[:-1,0:1] -1 instead of just all since the last row is rowindeces
    annObs = pd.DataFrame(index=df.T.iloc[:-1,0:1].index, data={'CellBarcode':df.T.iloc[:-1,0:1].index})

    #annObs = pd.DataFrame(index=completeTable.T.iloc[:,0:1].index, data={'CellBarcode' : .T.iloc[:,0:1].index,'Sample' : completeTable.T['0_Sample'].tolist()})

    annObs
    annObs.index
    annVar = pd.DataFrame(index=df.iloc[:,0:1].index, data=df.iloc[:,0:1].index, columns=['Gene'])
    annVar.index
    #adata = ad.AnnData(annMatrix, obs=annObs)

    #the levels are ogranised differently between annmatrix and annobs
    adata_TX = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
    adata_TX.obs_names_make_unique(join="-")
    #adata_TX.write(filename="/media/data/AtteR/scifi-analysis/Python-scifi-analysis/plots_no_groups/All.h5ad", compression=None, compression_opts=None, force_dense=None, as_dense=())
    return(adata_TX)

Pool1_final_RT_Hpos_full
Pool1_final_RT_Hneg_full
Pool2_RT_Hpos_full
Pool2_final_RT_Hneg_full



Pool1_RT_Hpos_adata = df_to_anndata(Pool1_final_RT_Hpos_full)
Pool1_RT_Hpos_adata
Pool1_RT_Hneg_adata = df_to_anndata(Pool1_final_RT_Hneg_full)
Pool1_RT_Hneg_adata
Pool2_RT_Hpos_adata = df_to_anndata(Pool2_RT_Hpos_full)
Pool2_RT_Hpos_adata
Pool2_RT_Hneg_adata = df_to_anndata(Pool2_final_RT_Hneg_full)
Pool2_RT_Hneg_adata

###################################################################
###################################################################
def preprocess(adata):
    #value counts returns the count of unique values
    print(adata.obs['CellBarcode'].value_counts())

    p = "/media/data/AtteR/scifi-analysis/scifi6/scanpy_analysis/plots/violinplot_P1Hp.pdf"

    adata.obs['CellBarcode']
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],multi_panel=True)

    sc.pl.scatter(P1Hp, "total_counts", "n_genes_by_counts")

    #knee plot
    #@title Threshold cells according to knee plot { run: "auto", vertical-output: true }
    #expected_num_cells =  50#@param {type:"integer"}
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, range(len(knee)), linewidth=5, color="g")
    #ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    #ax.axhline(y=expected_num_cells, linewidth=3, color="k")
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    plt.grid(True, which="both")
    plt.show()

    adata.X
    ##########
    gene_counts = np.sum(adata.X, axis=1)
    gene_counts.shape

    all_genes = adata.obs['total_counts']
    all_genes.shape
    print("Minimum number of transcripts per cell:", np.min(all_genes), 
        "\n Median number of transcripts per cell:", np.median(all_genes),
        "\n Maximum number of transcripts per cell:", np.max(all_genes))
    plt.hist(all_genes, bins=100)
    plt.savefig('SciFiAll_Histogram_2200.pdf')
    plt.show()
    return(adata)


P1Hp = Pool1_RT_Hpos_adata.copy()
P1Hn = Pool1_RT_Hneg_adata.copy()
P2Hp = Pool2_RT_Hpos_adata.copy()
P2Hn = Pool2_RT_Hneg_adata.copy()

P1Hp_A = preprocess(P1Hp)
P1Hn_A = preprocess(P1Hn)
P2Hp_A = preprocess(P2Hp)
P2Hn_A = preprocess(P2Hn)

#################

# Removes cells with less than X cells

adata = P1Hn_A[np.asarray(P1Hn_A.X.sum(axis=1)).reshape(-1) > 600]   #60,  640, 640

print("Cells prior to filtering %d"%P1Hn_A.n_obs)

print("Remaining cells %d"%adata.n_obs)

# Removes genes with 0 umi counts
adata = adata[:, np.asarray(adata.X.sum(axis=0)).reshape(-1) > 0]
adata

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],multi_panel=True)
# mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
# ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pl.highest_expr_genes(adata, n_top=20)
# filter for percent mito

adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# filter for percent ribo > 0.05
adata = adata[adata.obs['pct_counts_ribo'] > 5, :]

#sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],multi_panel=True)
# Create a mask to filter out cells with more than 6500 genes, less than 200 genes or less than 0.2 mitochondrial umi counts
mask = np.logical_or((adata.obs.n_genes < 1000).values, (adata.obs.n_genes > 0).values, (adata.obs.percent_mito < 0.2).values)
#Normalization
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4) # scale each cell to a common library size
sc.pp.log1p(adata) # log(expression + 1)
adata.raw = adata # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    # subset=True, # to automatically subset to the 4000 genes
    layer="counts",
    flavor="seurat_v3"
)
# import subprocess

# def install(name):
#     subprocess.call(['pip3', 'install', name])
# install("scikit-misc")
adata = adata[:, adata.var.highly_variable].copy()


#count number of umis per barcode
umi.per.barcode <- colSums(full_cm)
x <- sort(umi.per.barcode, decreasing = TRUE)
plot(x, log="xy",type="l", xlab="Barcodes", ylab="UMI counts")


P1Hp.layers["counts"] = P1Hp.X.copy() # preserve counts
sc.pp.normalize_total(P1Hp, target_sum=1e4) # scale each cell to a common library size
sc.pp.log1p(P1Hp) # log(expression + 1)
P1Hp.raw = P1Hp # freeze the state in `.raw`

sc.pp.highly_variable_genes(P1Hp, n_top_genes=4000, layer="counts", flavor="seurat_v3")
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
