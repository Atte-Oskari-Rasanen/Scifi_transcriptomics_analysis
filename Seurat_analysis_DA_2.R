library("Seurat")
library("dplyr")
library("patchwork")

library("Matrix")
library("readr")

setwd("/media/data/AtteR/scifi-analysis/R-scifi-analysis")
library ("berryFunctions")
library("fs")
library("plyr")

#Friday 04/02/22
#1. Process scifi5 only by using the marker genes, check if any matches and check their specificity for the certain clusters etc. redo clustering.Use lenient preprocessing so e.g. that each cell should contain min. 2400 genes. Get the 
# top marker genes from the list of the known markers only. 

#2. Explore ambient RNA removal methods

#3. Redo this but with all data, use more lenient filtering. Or process this data
# as separate from scifi5


#For generating the knee plot to know at which point the quality of cells 
#end up being poor. Afterwards the data is imported using Read10x to get 
#everyting as S4 objects instead of lists since if you combine the lists,
#you end up losing the gene names. 



#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#pbmc.data <- read.csv("/home/atte/Documents/Thesis/Scifi_library_2_S2_SciFi5_oDT/count_matrix.csv")
cm5_filt <- Read10X(data.dir = "/media/data/AtteR/scifi-analysis/outputs_starsolo/Scifi_library_2_S2_SciFi5_oDT/GeneFull/filtered")
cm5_filt <- Read10X(data.dir="/media/data/AtteR/scifi-analysis/outputs_starsolo/Scifi5/Scifi_library_2_S2_SciFi5_WP/GeneFull/filtered")

cm5_raw <- Read10X(data.dir = "/media/data/AtteR/scifi-analysis/outputs_starsolo/Scifi_library_2_S2_SciFi5_oDT/GeneFull/raw")

counts <- read_10x("/media/data/AtteR/scifi-analysis/outputs_starsolo/Scifi_library_2_S2_SciFi5_oDT/GeneFull/filtered") # Read 10X data into sparse matrix
obj_f

obj_f = CreateSeuratObject(counts = cm5_filt, min.cells = 3, min.features  = 200, project = "Scifi5", assay = "RNA")
table(obj_f$orig.ident)


library("diem")
# DIEM steps
sce <- set_debris_test_set(obj_f)
sce <- filter_genes(obj_f)
sce <- get_pcs(obj_f)
sce <- init(obj_f)
sce <- run_em(obj_f)
sce <- assign_clusters(obj_f)
sce <- estimate_dbr_score(obj_f)

# Evaluate debris scores
sm <- summarize_clusters(obj_f)
plot_clust(sce, feat_x = "n_genes", feat_y = "score.debris", 
           log_x = TRUE, log_y = FALSE)
plot_clust(sce, feat_x = "pct.mt", feat_y = "score.debris", 
           log_x = TRUE, log_y = FALSE)

# Call targets using debris score for single-nucleus data
sce <- call_targets(obj_f, thresh_score = 0.5)

# Call targets by removing droplets in debris cluster(s) for single-cell data
sce <- call_targets(obj_f, clusters = "debris", thresh = NULL)

seur <- convert_to_seurat(obj_f)
pbmc5 <- CreateSeuratObject(counts = cm5, min.cells = 3, min.features  = 200, project = "Scifi5", assay = "RNA")



#Only include the rows that contain the marker genes
DA_markers <- c('ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4', 'GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B', 'HOMER2','AJAP1') #genes from sncg are from paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6362095/#MOESM1   https://www.frontiersin.org/articles/10.3389/fcell.2020.00443/full 

DA_markers <- c('ABCC8','ADCYAP1','ADRA1B','AIF1','ALDH1A1','ALDOC','ANXA1','BNC2','C3','CADPS2','CALB1','CALB2','CALCA','CBLN1','CBLN4','CCK','CD9','CDH8','CHRNA5','CHST8','CLIP2','CNR1','CRHBP','CPLX1','DDC','DRD2','EN1','EN2','ENO3','EPHA4','ERLEC1','ETV1','EXT2','FGF1','FJX1','FOXA1','FOXA2','TH','SLC6A3','KCNJ6','CALB1', 'NR4A2', 'DRD2','GFAP', 'SNCG', 'VIP', 'NXPH4', 'GAD2', 'DAT', 'SOX6','VGLUT2', 'OTX2', 'CORIN1', 'WNT1', 'LMX1A', 'PBX1', 'PITX3', 'GRIN2B', 'HOMER2','AJAP1', 'EPHA4', 'CHRNA5', 'NRIP3', 'KCNS3', 'CPLX1', 'NDNF', 'KIFC3', 'CALB1', 'CHST8', 'IGFBP2', 'LAMA5', 'ANXA1', 'RBP4', 'ALDH1A7', 'ADCYAP1', 'LHFPL2', 'CBLN4', 'LPL', 'NHNIH2, OTX1, SYN2, CLBN1', 'GPX3', 'FJX1', 'FOXA2', 'EN2', 'NTF3', 'GFRA2', 'LIX1', 'PTPN5', 'FGF1', 'NOSTRIN', 'SERPINE2', 'KCNIP3', 'GRIK1', 'LYPD1', 'POU3F1', 'CD9', 'NEUROD6', 'GRP', 'TCF12', 'CALCA', 'GPR83', 'NPHP1', 'CHTF8', 'SLC32A1', 'CTXN3', 'ETV1', 'LMX1A')
#pbmc.big <- CreateSeuratObject(counts = full_cm, min.cells = 3, min.features  = 2400, project = "Scifi", assay = "RNA")

pbmc.big = obj_f

library("stringr")
DA_markers <- tolower(DA_markers)
DA_markers <- str_to_title(DA_markers)

#pbmc.DA = pbmc.big
#To get the dopaminergic marker genes only
#pbmc.DA <- subset(x=pbmc.big_f, features = DA_markers) #this was an earlier ver

pbmc.DA <- subset(x=pbmc.big, features = DA_markers)

length(rownames(pbmc.big))
length(rownames(pbmc.DA))

table(pbmc.DA_raw$orig.ident)

pbmc.DA <- subset(x=obj_f, features = DA_markers)

table(pbmc.DA_f$orig.ident)
#View(pbmc.big@meta.data)
# Add number of genes per UMI for each cell to metadata, i.e. number of genes found inside a cell divided # by total number of molecules inside the cell?

#number of genes detected per UMI
pbmc.DA$log10GenesPerUMI <- log10(pbmc.big$nFeature_RNA) / log10(pbmc.big$nCount_RNA)

# Create metadata dataframe
metadata <- pbmc.DA@meta.data
head(metadata)

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
head(metadata)

library("stringr")
# Create sample column
#metadata$sample <- NA
#metadata$sample[which(str_detect(metadata$cells, "^Scifi5_"))] <- "Scifi5"
#metadata$sample[which(str_detect(metadata$cells, "^Scifi6_"))] <- "Scifi6"
#metadata$sample[which(str_detect(metadata$cells, "^Scifi7_"))] <- "Scifi7"
#metadata$sample[which(str_detect(metadata$cells, "^Scifi8_"))] <- "Scifi8"

# Visualize QC metrics as a violin plot
VlnPlot(pbmc.DA, features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI"), ncol = 2)


# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.
#par(mfrow = c(1, 2))
#nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules #detected within a cell
plot2 <- FeatureScatter(pbmc.DA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1
plot2
#more genes åer cell compared total number of molecules per cell... does not make sense
#https://www.biostars.org/p/9493884/ 

#transcript count per cell
hist(pbmc.DA$nFeature_RNA, breaks=100)
hist(pbmc.DA$nCount_RNA, breaks=100)


counts <- GetAssayData(pbmc.DA, slot="counts", assay="RNA")   
genes.percent.expression <- rowMeans(counts>0 )  
genes.percent.expression
#Filter out cells with unique gene counts
# We filter out cells that have unique gene counts (nFeature_RNA) over 2,500 or less than
# 200 Note that > and < are used to define a'gate'.  
#-Inf and Inf should be used if you don't want a lower or upper threshold.
#pbmc.big_f <- subset(x = pbmc.big, subset = nFeature_RNA > 200 & nFeature_RNA < 5800)

#with cm5 only



#normalise data
pbmc.DA_f <- NormalizeData(pbmc.DA, normalization.method = "LogNormalize", scale.factor = 10000)

#Find variable genes, ignore housekeeping genes

#Approach 1
#Housekeeping genes that are similarly expressed in all cell populations are not useful for the purpose #of identifying these populations. Thus, it is often useful to select a subset of genes that display #higher than average variability among cells to be used for dimensionality reduction and clustering of #cells, as this will greatly speed-up the computations.
#base the nfeatures nbased on the scatter plot outliers, i.e. here we return 5500 genes per cell/dataset
pbmc.DA_f <- FindVariableFeatures(object = pbmc.DA_f, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 49)  #earlier 49
head(rownames(pbmc.DA_f))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc.DA_f), 10)
top10


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.DA_f)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
#view the output
head(x = HVFInfo(object = pbmc.big_f))




#data scaling and removing unwanted sources of variation like tehc. noise and batch effects

# scale data via linear transformation prior to pca
all.genes <- rownames(pbmc.DA_f)


pbmc.DA_f <- ScaleData(pbmc.DA_f, features = all.genes)


#The idea of this soupx method is that genes that are highly expressed in the soup and are marker genes for some population can be used to estimate the background contamination

#but we dont know the clusters markers yet....

head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)
#Scree plot prior to pca
#PCA  #had to increase nfeatures in findvarfeatures and when generating the seurat object...
#only keep the top 30 PCs 

#add approx false since earlier it gave a warning message...
pbmc.DA_f <- RunPCA(pbmc.DA_f, features = VariableFeatures(object = pbmc.DA_f))

#pbmc.DA_f$seurat_clusters
#sc = autoEstCont(sc)
# Examine and visualize PCA results a few different ways
print(pbmc.DA_f[["pca"]], dims = 1:3, nfeatures = 20)
DimPlot(pbmc.DA_f, reduction = "pca")

DimHeatmap(pbmc.DA_f, dims = 1, cells = 500, balanced = TRUE)

#Trying to define the number of dims via jacstraw and elbowplot
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values #for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment #of #features with low p-values (solid curve above the dashed line).
JackStrawPlot(pbmc, dims = 1:20)  # maybe draw the line at 10 pcs?

ElbowPlot(pbmc.DA_f, ndims = 23, reduction = 'pca') #between 9 and 12 pcs?

#so maybe at around 4?

#CLUSTER PCA
library(reticulate)
#py_install("umap-learn")
#reticulate::install_miniconda() 
#reticulate::py_install(packages = 'umap-learn')

####Harmony step (optional)
#library("harmony")
#pbmc <- RunHarmony(pbmc, "orig.ident")
#####

pbmc.DA_f <- FindNeighbors(pbmc.DA_f, dims = 1:22)
pbmc.DA_f <- FindClusters(pbmc.DA_f, resolution = 1.2)
pbmc.DA_f <- RunUMAP(pbmc.DA_f, reduction = "pca", dims = 1:22)
DimPlot(pbmc.DA_f, reduction = "umap", split.by = "seurat_clusters")
DimPlot(pbmc.DA_f, reduction = "umap")

source("chooseR.R")
library("clustree")

# plot the tree of clusters
# K1 to K5 correspond to different clustering resolutions
clustree(pbmc.DA_f$seurat_clusters, prefix = "RNA")

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc.DA_f), 5)

saveRDS(pbmc.DA_f, file = "/media/data/AtteR/scifi-analysis/R-scifi-analysis/DA_scifi5.rds")

##Finding markers
#Let’s make a “SoupChannel”, the object needed to run SoupX. Detailed info is available in SoupX vignette (see link above).
library("SoupX")
soup.channel  <- SoupChannel(cm5_raw, cm5_filt)
#SoupX requires clusters in order to define marker genes. You can either use CellRanger clustering (see SoupX vignette), or quickly cluster using Seurat. We’ll go for option 2 here. All the following steps will be addressed in more detail in the separate Seurat vignette.

pbmc.DA_f<- SCTransform(pbmc.DA_f, verbose = F)
pbmc.DA_f<- RunPCA(pbmc.DA_f, verbose = F)
pbmc.DA_f<- RunUMAP(pbmc.DA_f, dims = 1:22, verbose = F)
pbmc.DA_f<- FindNeighbors(pbmc.DA_f, dims = 1:22, verbose = F)
pbmc.DA_f<- FindClusters(pbmc.DA_f,resolution = 1.2, verbose = T)

meta<- pbmc.DA_f@meta.data
umap<- pbmc.DA_f@reductions$umap@cell.embeddings
soup.channel<- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel<- setDR(soup.channel, umap)
head(meta)

#With defined clusters, run the main SoupX function, calculating ambient RNA profile.
#quickMarkers(cm5_filt, pbmc.DA_f$seurat_clusters, N = 10, FDR = 0.01, expressCut = 0.9)
#an issue 
soup.channel  <- autoEstCont(soup.channel, tfidfMin=0.9)

# find markers for every cluster compared to all remaining cells, report
# only the positive ones

# findallmarkers compares the cluster at hand to all other clustes so may not
# be that informative. need to compare individual ones. so if we would like to determine the genes that are differentially expressed between these specific clusters, we can use the FindMarkers() function


library("scSorter")
markers <- FindAllMarkers(object = pbmc.DA_f, only.pos = TRUE, min.pct = 0.2,min.diff.pct=0.3, thresh.use = 0.5)
#Q
#Any point in developing our own classifier based on dopaminergic cell subtypesin case we dont have a pre-existing dataset that has this already?

#avg_logFC: average log2 fold change. Positive values indicate that the gene is # more highly expressed in the cluster.
#pct.1: The percentage of cells where the gene is detected in the cluster
#pct.2: The percentage of cells where the gene is detected on average in the other clusters
markers <- markers[ markers$p_val_adj < 0.01, ]
head(markers)
View(markers)
#View(top10_markers)
write.table(markers, file="/media/data/AtteR/scifi-analysis/R-scifi-analysis/markers/DA_only_markers_res02_pcs3.csv")

#should I look at the marker that it most significant even if its pct values wouldnt be very different?
markers_0vs1 <- FindMarkers(object = pbmc.DA_f, ident.1 = 0, ident.2 = 1)
markers_0vs1
#0: kifc3, cck, ndnf, kcnip3, ntf3, gfra2, chrna5, foxa1
# chrna5, ndnf, kcns3, en1
markers_0vs2 <- FindMarkers(object = pbmc.DA_f, ident.1 = 0, ident.2 = 2)
#2: Anxa1, Grp, Foxa2, Chst8

markers_0vs3 <- FindMarkers(object = pbmc.DA_f, ident.1 = 0, ident.2 = 3)

#3: Gfap, Scl32a1

markers_0vs4 <- FindMarkers(object = pbmc.DA_f, ident.1 = 0, ident.2 = 4)
#4 ctxn3, c3, Aldoc
markers_0vs5 <- FindMarkers(object = pbmc.DA_f, ident.1 = 0, ident.2 = 5)
#5: Pou3f1, Gfap

features_r04_pcs22 <- c("Kifc3", "Cck", "Ndnf", "Kcnip3", "Ntf3", "Gfra2", "Chrna5", "Foxa1", "Anxa1", "Grp", "Foxa2", "Chst8", "Scl32a1", "Ctxn3", "C3", "Aldoc", "Pou3f1", "Gfap")
#for all
#VlnPlot(pbmc, features = c("RGD1566401", "Frmd8", "Tll1", "Frmd81", "Asic2", "Lrrc4c"))
#cm5 pcs 4

features_r02_pcs15 <- c("Calb1", "Gfap", "Foxa2")
VlnPlot(pbmc.DA_f, features = features_r08_pcs22)
#cm5 pcs3
VlnPlot(pbmc, features = c("Frmd8", "Cdc42ep2", "Frmd81", "Rbfox1", "Cntnap5a", "Asic2"))



markers.0 <- markers[ which(markers$cluster == 0), ]
markers.1 <- markers[ which(markers$cluster == 1), ]
markers.2 <- markers[ which(markers$cluster == 2), ]
markers.3 <- markers[ which(markers$cluster == 3), ]
markers.4 <- markers[ which(markers$cluster == 4), ]

#head(markers.1)


c1 <- c("Chrna5", "Ndnf", "Kcns3", "En1")
c2 <- c("Anxa1", "Grp", "Foxa2")
c3_7 <- c("Cbln4", "Anxa1", "Gfap", "Slc32a1", "Foxa2", "Ctxn3", "Pou3f1")

#2 ekaa 3, 2 tokaa 4, 2 seur 5, sit 6, 6
features_r08_pcs22 <- c(c1,c2,c3_7)
features_r04_pcs22_a <- c("Kifc3", "Cck", "Ndnf", "Kcnip3", "Ntf3", "Gfra2", "Chrna5", "Foxa1", "Anxa1", "Grp", "Foxa2")
features_r04_pcs22_b <- c("Chst8", "Scl32a1", "Ctxn3", "C3", "Aldoc", "Pou3f1", "Gfap")

features_r12_pcs22_c1 <- c("Aldh1a1", "Chrna5", "Lix1")
features_r12_pcs22_c2 <- c("Chst8", "Nr4a2")
features_r12_pcs22_c3 <- c("Cck", "Kifc3")
features_r12_pcs22_c5 <- c("Gfap", "Slc32a1")
features_r12_pcs22_c6 <- c("Calb1", "Ndnf", "Ndnf")

features_r12_pcs22_c7 <- c("Cbln4")
features_r12_pcs22_c8 <- c("Grp", "Foxa2")

features_r12_pcs22_c9_10_11 <- c("Ctxn3", "Anxa1", "Pou3f1")
features_r12_pcs22 <- c(features_r12_pcs22_c1,features_r12_pcs22_c2, features_r12_pcs22_c3, features_r12_pcs22_c5, features_r12_pcs22_c6, features_r12_pcs22_c7, features_r12_pcs22_c8, features_r12_pcs22_c9_10_11 )

FeaturePlot(pbmc.DA_f, features = features_r12_pcs22, reduction = "umap")

VlnPlot(pbmc.DA_f, features = features_r12_pcs22)

FeaturePlot(pbmc.DA_f, features = features_r08_pcs22, reduction = "umap")

DotPlot(pbmc.DA_f, features=features_r12_pcs22)


#Cluster renaming

pbmc.DA_f_rename <- RenameIdents(object = pbmc.DA_f, `0` = "0", `1` = "SNc", `2` = "DA-VTA2/3/4", `3`= "DTA-VTA3", `4`="4", `5`="DA-VTA4, SN", `6`="DA-VTA1", `7`="DA-VTA2/3", `8`= "DA-VTA2/1",`9`="DA-VTA4", `10`="DA-VTA2, SN", `11`= "DA-VTA1")

DimPlot(pbmc.DA_f_rename, reduction = "umap")

#Clip2 insig
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
features <- c("X", "Nr4a2", "Cplx1", "Clip2(insig)", "Cnr1")

features <- (c("RGD1566401", "Frmd8", "Tll1", "Frmd81", "Asic2", "Lrrc4c"))
RidgePlot(pbmc.DA_f, features = features_r02_pcs15, ncol = 2)

#5,6 appear to be quite similar so probs a result of overclustering?

#save the original clusters
pbmc0[["OriginalClusterNames"]] <- Idents(object = pbmc)

current.cluster.ids <- c(0,1,2,3,4,5)
new.cluster.ids <- c(0,1,2,3,45,45)
pbmc0$OriginalClusterNames <- plyr::mapvalues(pbmc0$OriginalClusterNames, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(pbmc0, reduction = "umap")

#ERROR
#pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Visualise top markers as heatmap
top.markers <- do.call(rbind, lapply(split(markers, markers$cluster), head))
#DoHeatmap(pbmc, genes.use = top.markers$gene, slim.col.label = TRUE, remove.key = TRUE)
markers.0 <- markers[ which(markers$cluster == 0), ]

VlnPlot(pbmc, features.plot = head(markers.0$gene), point.size.use=0.5)

## calculate cluster enriched genes for leiden clustering

# libraries
library(ggraph)
library(igraph)
library(tidyverse)

# create an edge list data frame giving the hierarchical structure of your individuals
d1 <- data.frame(from="origin", to=paste("group", seq(1,5), sep=""))
d2 <- data.frame(from=rep(d1$to, each=5), to=paste("subgroup", seq(1,25), sep="_"))
edges <- rbind(d1, d2)

# Create a graph object 
mygraph <- graph_from_data_frame( edges )

# Basic tree
ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal() +
  geom_node_point() +
  theme_void()
#different methods for var genes discovery also available
var_genes <- VariableFeatures(pbmc.big_f)
seurat_df <- GetAssayData(pbmc.big_f)[var_genes,]
