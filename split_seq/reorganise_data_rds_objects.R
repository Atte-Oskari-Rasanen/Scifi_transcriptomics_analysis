library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)



#install.packages('Seurat')


pbmc.data <- Read10X(data.dir = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


pbmc$RNA@data



rm(list = ls())

sample_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Merged/mouse_merged/hum_asyn1_sn_nogfp"
fig_path='/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Merged/plots/'

data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/rat/analysis_rat/all-well"
data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample2/mus/analysis_mus_s2/"

data_path="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/analysis_mus_s1/"


data_path_m3="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Merged/mouse_merged"
data_path_r3="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Merged/rat_merged"


data_path_m2="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample2/mus/analysis_mus_s2"
data_path_r2="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample2/rat/analysis_rat_s2"

data_path_r1="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/rat/analysis_humrat_s1/"
data_path_m1="/media/data/AtteR/projects/parsebio/ParsePipeline/sc_countmatrices/Sample1/mus/analysis_mus_s1"
data_paths=list(data_path_r1, data_path_r2, data_path_m2, data_path_r3, data_path_m3)
data_paths=list(data_path_r2, data_path_m2, data_path_r3, data_path_m3)
data_paths=list(data_path_m1)

data_paths=list(data_path_r3, data_path_m3)

#Mus sample 1, sample 2, merge indiv ones, then combine the animals of each group and rerun
#take into consideration DA vs rest and control ones

print ("list of directories including the main directory")
dir_list <- list.dirs(data_path_r2,recursive = FALSE) 
samples_trim <-Filter(function(x) !any(grepl("process", x)), dir_list)
samples_trim <-Filter(function(x) !any(grepl("plots", x)), samples_trim)




for (data_path in data_paths){
  dir_list <- list.dirs(data_path,recursive = FALSE) 
  dir_list
  ##############################
  #samples <- dir_list[grepl("mus",dir_list)] #6-gfp,7,8 -nogfp
  #da_samples <- samples[grepl("gfp",samples)]
  #rest_samples <- samples[grepl("nogfp",samples)]
  
  samples_trim <-Filter(function(x) !any(grepl("process", x)), dir_list)
  samples_trim <-Filter(function(x) !any(grepl("plots", x)), samples_trim)
  samples_trim <-Filter(function(x) !any(grepl("all-well", x)), samples_trim)
  
  #length(samples)
  length(samples_trim)
  
  for (sample_path in samples_trim){
    print(sample_path)
    fig_path<-paste(sample_path,"/plots", sep="")
    dir.create(file.path(sample_path, 'plots'), showWarnings = FALSE)
    
    ##############################
    ##############################
    
    mat <- readMM(paste0(sample_path, "/DGE_filtered/DGE.mtx"))
    print("matrix read")
    cell_meta <- read.delim(paste0(sample_path, "/DGE_filtered/cell_metadata.csv"),
                            stringsAsFactor = FALSE, sep = ",")
    genes <- read.delim(paste0(sample_path, "/DGE_filtered/all_genes.csv"),
                        stringsAsFactor = FALSE, sep = ",")
    
    cell_meta$bc_wells <- make.unique(cell_meta$bc_wells, sep = "_dup")
    
    rownames(cell_meta) <- cell_meta$bc_wells
    print("---")
    genes$gene_name <- toupper(genes$gene_name)
    
    genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")
    print("---")
    
    # Setting column and rownames to expression matrix
    colnames(mat) <- genes$gene_name
    rownames(mat) <- rownames(cell_meta)
    mat_t <- t(mat)
    
    df <- as.data.frame(as.matrix(mat))
    length(rownames(df))
    ######################################
    
    #get the different indeces
    copy_no=unique(df$AAV_ASYN)
    if (is.null(copy_no)){
      break      
    }
    #populate the vector of lists with the subsetted dataframes
    xy <- vector("list", length(copy_no))
    
    #subset the dataframe based on AAV_ASYN copy number
    i<-1
    for (copy in copy_no) {
      print(copy)
      xy[[i]]<- df %>% filter(df$AAV_ASYN == copy )
      i<-i+1
    }
    df_sorted <- do.call(rbind, xy)
    df_sorted[1:2,1:3]
    colnames(df_sorted)
    length(rownames(df_sorted))
    
    newdata <- df_sorted[c(length(rownames(df_sorted)),1)]
    head(newdata)
    
    #df_aav = subset(newdata, select = -c(Mmp14))
    df_aav <- newdata[c(0,-1)]
    
    head(df_aav)
    length(rownames(df_aav))
    # Remove empty rownames, if they exist
    #mat_t <- mat_t[(rownames(mat_t) != ""),]
    pbmc <- CreateSeuratObject(counts = t(df_sorted))
    # create a new assay to store AAV information
    aav_assay <- CreateAssayObject(counts = t(df_aav))
    pbmc[["AAV"]] <- aav_assay
    
    #SaveObject(pbmc, "seurat_obj_before_QC")
    saveRDS(pbmc, file = paste(sample_path,"/orig.RDS",sep=""))
    
  }
}

