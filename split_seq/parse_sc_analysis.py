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

pseudo code for merging the samples 
extract given df from sample 1, find the equivalent from sample 2, import as anndatas, make into dfs
and concatenate with pandas

'''
class data_tools():
    def __init__(self, data_path, sample, injection="asyn"):
        self.data_path=data_path
        self.sample=sample
        self.injection=injection
    def __str__(self):
        return("a class for performing simple processes regarding anndata formatting")

    def get_sample_paths(self):
        samples_list = [dir for dir in sorted(glob.glob(f'{self.data_path}*{self.sample}*')) if os.path.isdir(dir)] 
        if self.sample=="all-well":
            return(samples_list)
        else:
            if self.injection=="asyn":
                paths_sub=[]
                for p in samples_list:
                    if "asyn" in p:
                        paths_sub.append(p)
            elif self.injection=="bfp":
                paths_sub=[]
                for p in samples_list:
                    if "bfp" in p:
                        paths_sub.append(p)


        return(paths_sub)

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
        gene_cell_df=adata.to_df()
        return([gene_cell_df, adata.obs['gene_count'], adata.obs['tscp_count']])

    ########################################################################################################
    def intersect_genes(self, samples_list):
        all_genes=[]
        for sample_df in samples_list:
            all_genes.append(list(set(sample_df.columns)))     
        common_genes=list(set.intersection(*map(set,all_genes)))
        # check length
        if len(common_genes) > 0:
            return(common_genes) 
        else:
            return("no common elements")
    def df_to_anndata(self, df):
        #annMatrix = df_grouped.T.iloc[:,1:] #why was anndata earlier taken like this?
        #annMatrix = df #[6987 rows x 163 columns]
        annMatrix = df.iloc[1:,:] #[6987 rows x 163 columns]

        #annMatrix = df_grouped.T.iloc[:,1:] #why was anndata earlier taken like this?

        annObs = pd.DataFrame(index=df.T.iloc[:,0:1].index, data={'CellBarcode' : df.T.iloc[:,0:1].index, 'AAV_ASYN': df.iloc[0,:].tolist()})
        annVar = pd.DataFrame(index=df.iloc[1:,0].index, data=df.iloc[1:,0].index, columns=['Gene'])
        adata = ad.AnnData(X = annMatrix.T, obs = annObs, var = annVar)
        adata.obs_names_make_unique(join="-")
        return(adata)
    def grouping_aav(self, df,aav):
        aav=str(aav).upper()
        dfs_all=[]
        df_reorg = df.sort_values(by = aav)
        aav_copies=list(set(df_reorg[aav]))
        for aav_copy in aav_copies:
            aav_col=np.array(df_reorg[aav])
            #get the indices based on which you will subset the Df_reorg
            indeces = np.where(np.array(aav_col) == aav_copy)
            indeces = [item for t in indeces for item in t] #need to transform the tuple into a list
            df_aav=df_reorg.iloc[indeces[0]:indeces[-1],:]
            dfs_all.append(df_aav)
            print(f'extracted a sub array with aav_copy number {aav_copy} out of {len(aav_copies)}')
        df_sorted = pd.concat(dfs_all).fillna(0)
        #df_sorted = reduce(lambda df_left,df_right: pd.concat(df_left, df_right), dfs_all).fillna(0)
        return(df_sorted.T)
    def merge_all(self):
        samples_list=self.get_sample_paths()
        all_dfs_and_counts = list(map(self.retrieve_data, samples_list))
        all_dfs_sample=[]
        all_gene_counts=[]
        all_tcp_counts=[]
        for i in range(len(all_dfs_and_counts)):
            all_dfs_sample.append(all_dfs_and_counts[i][0])
            all_gene_counts.append(pd.Series.to_frame(all_dfs_and_counts[i][1]))
            all_tcp_counts.append(pd.Series.to_frame(all_dfs_and_counts[i][2]))
        gene_counts=pd.concat(all_gene_counts)
        tcp_counts=pd.concat(all_tcp_counts)
        common_genes = self.intersect_genes(all_dfs_sample)
        df_merged = reduce(lambda df_left,df_right: pd.merge(df_left, df_right,on=common_genes, how='outer'), all_dfs_sample).fillna(0)
        #remove x and y from the non-overlapping ones
        df_merged.drop(df_merged.filter(regex='_y$').columns, axis=1, inplace=True)
        df_merged.drop(df_merged.filter(regex='_x$').columns, axis=1, inplace=True)
        #adata = self.df_to_anndata(df_merged)
        med_g=str(gene_counts.median(0)).split("\n")[0].split("    ")[1]
        med_t=str(tcp_counts.median(0)).split("\n")[0].split("    ")[1]
        print('median transcript count per cell: ' + str(tcp_counts.median(0)))
        print('median gene count per cell: ' + str(gene_counts.median(0)))

        return(df_merged.T)

#DO YOU NEED TO GROUP BASED ON ASYN?
#try find the gene column for a_syn, if found, then need to group based on it...

data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/rat/analysis_rat/all-well"
data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample2/mus/analysis_mus_s2/"

data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/analysis_mus_s1/"
data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/rat/analysis_humrat_s1/"
dat_format=data_tools(data_path, 'rat')
paths=dat_format.get_sample_paths()
len(paths)
df=dat_format.merge_all() #warning:AnnData expects .var.index to contain strings, but got values like:[0,1,3]
df
aav=str('aav_asyn').upper()
dfs_all=[]
df_reorg = df.T.sort_values(by = aav)
aav_copies=list(set(df_reorg[aav]))
for aav_copy in aav_copies:
    aav_col=np.array(df_reorg[aav])
    #get the indices based on which you will subset the Df_reorg
    indeces = np.where(np.array(aav_col) == aav_copy)
    indeces = [item for t in indeces for item in t] #need to transform the tuple into a list
    df_aav=df_reorg.iloc[indeces[0]:indeces[-1],:]
    dfs_all.append(df_aav)
    print(f'extracted a sub array with aav_copy number {aav_copy} out of {len(aav_copies)}')
df_sorted = pd.concat(dfs_all).fillna(0)


df_grouped=dat_format.grouping_aav(df,'aav_asyn')
adata=dat_format.df_to_anndata(df_reorg.T)
df=df_sorted.T
annMatrix=df_sorted.iloc[:,1:]
annMatrix

#somethings going on with the indeces of the merged df when trying to convert into anndata....
###################################################
annObs = pd.DataFrame(index=df.T.iloc[:,0:1].index, data={'CellBarcode' : df.T.iloc[:,0:1].index, 'AAV_ASYN': df.iloc[0,:].tolist()})
annVar = pd.DataFrame(index=df.iloc[1:,0].index, data=df.iloc[1:,0].index, columns=['Gene'])
adata = ad.AnnData(X = annMatrix, obs = annObs, var = annVar)
adata.obs_names_make_unique(join="-")
###################################################

#adata=dat_format.retrieve_data(paths[1])
#adata=sc.AnnData(adata[0])

#make an option to not only merge the sample mus ones but also to be able to merge sample1 and 2

'''
########################################
#save a holistic csv file where you get the median N(tcp)/cell and median N(gene)/cell for all samples
#for the mouse: merged and then each individual one
#calculate
med_g=str(gene_counts.median(0)).split("\n")[0].split("    ")[1]
med_t=str(tcp_counts.median(0)).split("\n")[0].split("    ")[1]

print("Minimum number of transcripts per cell:", np.min(all_genes), 
    "\n Median number of transcripts per cell:", np.median(all_genes),
    "\n Maximum number of transcripts per cell:", np.max(all_genes))

prel_info = [[med_g,med_t]]
df = pd.DataFrame(prel_info, columns=['median transcript count per cell','median gene count per cell', 'Minimum number of transcripts per cell', 'Maximum number of transcripts per cell'])
results_file = '/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/anndata_objects/'  # the file that will store the analysis results
df.to_csv(results_file+"tcp_gene_medians_cell.csv", encoding='utf-8', index=False)
########################################
'''

#median transcript count/cell
#median gene count/cell-so calculate counts for each cell, i.e. row, make a final sum column and then take the median
'''
          Gmpr  Rtca  Tbrg4  Ndufa9  Nsun5
05_01_96   0.0   1.0    0.0     0.0    0.0
05_03_09   0.0   0.0    0.0     0.0    0.0
05_03_58   0.0   0.0    0.0     0.0    0.0
05_03_96   0.0   1.0    0.0     1.0    0.0
05_04_09   2.0   0.0    0.0     0.0    0.0

1. does the data frame even include transcript counts? 
2.iterate over each row(cell), get the number of gene counts, make a list based on that. 

maybe gene count is actually the numbner of genes that are found in the cell and transcript
count is the number of transcripts mapped to the specific gene?

adata.uns["gene_count"] = 
'''

#DONT forget the control bfp sample!
genes=["AAV_CTE", "AAV_ASYN", "AAV_H2BGFP", "AAV_BFP"]
adata=dat_format.retrieve_data(paths[5])

#no asyn found in rat paths[1]- 'humrat_asyn1_sn_nogfp' 
df=adata[0] #[166 rows x 6987 columns]
df_grouped=dat_format.grouping_aav(df,'aav_asyn')

df=df_grouped
annMatrix = df_grouped.iloc[1:,:] #[6987 rows x 163 columns]
#annMatrix = df_grouped.T.iloc[:,1:] #why was anndata earlier taken like this?
annMatrix=df.iloc[1:,:]
annMatrix

###################################################
annObs = pd.DataFrame(index=df.T.iloc[:,0:1].index, data={'CellBarcode' : df.T.iloc[:,0:1].index, 'AAV_ASYN': df.iloc[0,:].tolist()})
annVar = pd.DataFrame(index=df.iloc[1:,0].index, data=df.iloc[1:,0].index, columns=['Gene'])
adata = ad.AnnData(X = annMatrix.T, obs = annObs, var = annVar)
adata.obs_names_make_unique(join="-")
###################################################

adata=dat_format.df_to_anndata(df_grouped)
adata.var['Gene']
#based on the grouped values we add an extra row to the top of the genes which shows the number of 
#asyn molecules in each cell




#annMatrix['redundant']=[0.0]*len(df_grouped.columns)
#annVar = pd.DataFrame(index=df_grouped.T.iloc[1:,0:1].index, data=df_grouped.T.iloc[1:,0:1].index, columns=['Gene'])

#final_df_T_g    # had to add final_df_T_WORKS.T.iloc[:-1,0:1] -1 instead of just all since the last row is rowindeces
#annObs = pd.DataFrame(index=final_df_T_WORKS.T.iloc[:-1,0:1].index, data={'CellBarcode':final_df_T_WORKS.T.iloc[:-1,0:1].index, 'N_asyn':final_df_T_WORKS.loc['ASYN_COPIES']})




#should i make separate anndata objects for all?

#now we have dfs for each individual group. Should we have instead used the DF_All_T_reorg and just use the copy number groups as equivalents of "samples" used prior?
DA_markers = ['TH', 'DAT', 'A1BG', 'ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4','GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B','HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8','IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2', 'OTX1', 'SYN2', 'CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83','NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A']
##############################

#results_file = '/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/anndata_objects/'  # the file that will store the analysis results
#adata.write(results_file+'mus_rat_sn.h5ad', compression="gzip")

fig_path='/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/rat/plots/'

##############################
sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')
sc.settings.figdir = fig_path
##############################

#save file format: [Sample]_[Animal]_[AAV+copy_number]_[Brain_area]_[GFP_state]_[Figure_type].png

sample='1'
animal='r6'
aav='asyn'
##############
#Quality control
##############
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata.obs['n_genes_by_counts']
adata.obs['total_counts']
adata.var['Gene']

# Scanpy will prepend the string in the save argument with "violin"
# and save it to our figure directory defined in the first step.
'''
the number of genes expressed in the count matrix (so how often a certain gene is expressed in the matrix?)
the total counts per cell
the percentage of counts in mitochondrial genes
'''

sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, save=f'{sample}_{animal}_{aav}_ngenes_counts.png') #gene counts is 1800
#plt.savefig(fig_path+, bbox_inches="tight")

sc.pl.violin(adata, ['total_counts'], jitter=0.4, save=f'{sample}_{animal}_{aav}_total_counts.png') #total counts is 4000...

#Filtering
# Filter the data
adata = adata[adata.obs.n_genes_by_counts < 7000,:] #merged:800
adata = adata[adata.obs.total_counts < 38000,:] #merged:1500
#adata = adata[adata.obs.pct_counts_mt < 15,:]
adata.shape # Checking number of cells remaining
sc.pl.highest_expr_genes(adata, n_top=20, )

#Visualise
#sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save=f'{sample}_{animal}_{aav}_genes_x_total_counts.png')

all_genes = adata.obs['total_counts']
all_genes.shape
print("Minimum number of transcripts per cell:", np.min(all_genes), 
    "\n Median number of transcripts per cell:", np.median(all_genes),
    "\n Maximum number of transcripts per cell:", np.max(all_genes))

#Normalise
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
sc.pl.highly_variable_genes(adata) # scanpy generates the filename automatically

# Save raw expression values before variable gene subset
adata.raw = adata

adata_orig=adata.copy()

adata = adata[:, adata.var.highly_variable]

#batch effect correction:
#a simple batch effect correction
sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)


'''
ComBat assumes that the biol. state of each cell is already known and thus it separates 
biol.effects from batch effects via a linear model --- BUT with scRNA data the cell type identity
of individual cells may not always be known...

mnmCorrect - applies mutual nearest neighbours between celss in different batches --> this way
identifies the common biol.condititions across batches post hoc
-->removes batch effects from the gene expression matrix using PCA
'''



#ComBat 
#sc.pp.combat(adata, ['AAV_ASYN'])

# mnn_correct
#sc.external.pp.mnn_correct(adata, ['total_counts'])

#PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50) # scanpy generates the filename automatically

#UMAP+LEIDEN
'''
sc.pp.neighbors=Compute a neighborhood graph of observations [McInnes18].
The neighbor search efficiency of this heavily relies on UMAP [McInnes18], which also provides a method for estimating connectivities of data points - the connectivity of the manifold (method=='umap'). If method=='gauss', connectivities are computed according to [Coifman05], in the adaption of [Haghverdi16].
'''
#adata.obs['AAV_ASYN']
#adata.var_names
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=7)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.2)
#sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save=f'{sample}_{animal}_{aav}_umap.png')
sc.pl.umap(adata, color=['leiden', 'AAV_ASYN'], frameon=False, save=f'{sample}_{animal}_{aav}_umap.png')
DA_markers = ['TH', 'DAT', 'A1BG', 'ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4','GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B','HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8','IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2', 'OTX1', 'SYN2', 'CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83','NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A']


#the given DA markers
import seaborn as sns
#sc.pl.umap(adata, color=DA_markers, color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True), ncols=3, vmax=25, frameon=False)
#Finding the cluster markers
#we will for now use the simple t-test but others are also found here:https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html#Finding-marker-genes
sc.tl.rank_genes_groups(adata.var['Gene'], 'leiden', method='t-test')


# The head function returns the top n genes per cluster
top_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
print(top_markers)
results_file='/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/rat/anndata_objects/'
adata.write(results_file + 'Sample1_Rat_R5.h5ad')

adata=ad.read_h5ad(results_file + 'Sample1_Rat_all.h5ad')
