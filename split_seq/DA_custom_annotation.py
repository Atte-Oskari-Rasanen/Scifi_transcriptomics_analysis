import os
import pandas as pd
loom_dir="/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref"

pbmc_neur_path = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/rds_objects"
os.environ["HDF5_USE_FILE_LOCKING"]

#ref = pd.read_table(loom_dir+"/PanglaoDB_markers_27_Mar_2020.tsv")
ref = pd.read_csv(pbmc_neur_path+"/all_pbmc.csv")
ref.shape
ref.index
len(ref.columns)
ref.iloc[0,0]
ref.columns[0]

#gene_names = list(list(zip(*ref.iloc[:,0]))[0])

gene_names=[]
counts = [] #a list of counts per gene for all cells
for g in ref.iloc[:,0]:
    #print(g.split(" ")[0])
    gene_names.append(g.split(" ")[0])
    counts.append(g.split(" ")[0:])

ref.iloc[:,0]

len(counts)
len(gene_names)
len(counts[0])
len(ref.iloc[:,0])

cell_ids = str(ref.columns).split(" ")[:-1]
#cell_ids = str(ref.columns).split(" ")
len(cell_ids)
cell_ids[0]=str(cell_ids[0].split("(")[1].replace("['", ""))
cell_ids[-1]=cell_ids[0].split("\\")[0].replace("Index(['","")
cell_ids[-1]

def prune(s):
    return(s.strip('\"')) 
cell_ids = list(map(prune, cell_ids))

len(cell_ids)
len(gene_names)

df = pd.DataFrame(counts)
len(df.columns) #cells
len(df.index) #genes


team.columns

#reconstructing the count matrix. now we name the cells based on indeces of the annotated df
#we label the cells with the new names... how do we transfer this metadata to the seurat object?
df = pd.DataFrame(counts, index =gene_names,
            columns =cell_ids)

df_da_ann.iloc[:,0]


# if the certain cell type is found in more than one, e.g. DA_A9 + DA_SNc
# maybe create these double identities first before annotating cells?


set(df_da_ann["DA_A9"])
set(df_da_ann["DA_SNc"])

for cell_type in df_da_ann.columns:
    print(cell_type)
    set(df_da_ann[cell_type])
    for cell_i in df_da_ann[cell_type]:
        if cell_i!='NaN':


'''
Idents(object = pbmc, cells = da_a9_i) <- "orig.ident"
Idents(object = pbmc, cells = da_snc_i) <- "orig.ident"
Idents(object = pbmc, cells = da_vta1_i) <- "orig.ident"
Idents(object = pbmc, cells = da_vta2_i) <- "orig.ident"
Idents(object = pbmc, cells = da_vta3_i) <- "orig.ident"
Idents(object = pbmc, cells = da_vta4_i) <- "orig.ident"
'''


df.iloc[:,1]


#########
#rename given columns based on the indeces of annotated DA cells, e.g. DA_A9_1, DA_A9_2 ...
#########

#Filter df, initally get only the ones htat contain Th > 0 
for col_i, cell in enumerate(df_da.columns):
    found_genes_cell=[]
    for row_i in range(len(df_da.index)-1):
        if df_da.index[row_i] ="Th":
            if float(df_da.iloc[row_i,col_i])>0:




#rename columns

da_genes =["TH", "DDC", "GCH1", "SCL18A2", "SLC6A3","EPHA4", "CHRNA5", "NRIP3", 
                                    "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
                                    "LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3",
          "FJX1", "FOXA2", "EN2", "NTF3", "GFRA2", "KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", 
        "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA", "NTF3", "GFRA2", "VIP", "CCK",
         "SYN2", "CBLN1", "GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]

def reformat(name):
    return(name.lower().capitalize())
da_genes = list(map(reformat, da_genes))
da_genes

da_genes_found = list(set(da_genes) & set(list(df.index)))
len(da_genes_found)
df_da=df.loc[da_genes_found,:]
df_da


DA_A9 = ['Th', 'Ddc', 'Gch1', 'Slc18A2']
DA_A9 = ['Th', 'Slc6a3', "EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
"LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3"]
DA_SNc = ["EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
"LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3"]
DA_VTA1 = ["EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF",
             "FJX1", "FOXA2", "EN2", "NTF3", "GFRA2"]
DA_VTA2 = ["KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", "ALDH11A", "ALDH1A7", 
            "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA"]
DA_VTA3 = ["KLFC3", "CALB1", "CHST8", "CBLN4", "LPL", "NHLH2", "OTX1", "SYN2",
             "CBLN1", "GPX3", "FJX1", "FOXA2", "EN2","NTF3", "GFRA2", "VIP", "CCK"]
DA_VTA4 = ["KLFC3", "CALB1", "CHST8", "SYN2", "CBLN1","GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]

cell_types_da={'DA_A9':['Th', 'Slc6a3', "EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
"LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3"], 'DA_VTA1' : ["EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF",
"FJX1", "FOXA2", "EN2", "NTF3", "GFRA2"],
'DA_VTA2' : ["KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", "ALDH11A", "ALDH1A7", 
            "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA"],
'DA_VTA3' : ["KLFC3", "CALB1", "CHST8", "CBLN4", "LPL", "NHLH2", "OTX1", "SYN2",
             "CBLN1", "GPX3", "FJX1", "FOXA2", "EN2","NTF3", "GFRA2", "VIP", "CCK"],
'DA_VTA4' : ["KLFC3", "CALB1", "CHST8", "SYN2", "CBLN1","GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]}
#######################
base_markers_da={'DA_A9':['TH', 'DDC', 'GCH1', 'VMAT', 'SLC6A3']}


#for c in cell_types_da.keys():
#    da_genes_found = list(set(da_genes) & set(cell_types_da[c]))
#    cell_types_da[c]=da_genes_found
for i in cell_types_da.keys():
    cell_types_da[i]=da_genes = list(map(reformat, cell_types_da[i]))

cell_types_da["DA_A9"]

df.loc[DA_A9]
df_da.loc['Th']


gfp_cols = [i for i in df_da.columns if "GFP+" in i]
gfp_cols

df_da_sub=df_da[gfp_cols]
df_da_sub
#for gfp in df_da.columns:
#    if "GFP+" in gfp:


da_type_coord = []
da_types_cell={'DA_A9': [], 'DA_VTA1' : [],
'DA_VTA2' : [],
'DA_VTA4' : []}
da_type_coord=[]


def determine_cell_type(cell_type_genes, df_da):
    da_type_coord = []
    da_types_cell={}
    print(cell_types_da)
    for col_i, cell in enumerate(df_da.columns):
        found_genes_cell=[]
        for row_i in range(len(df_da.index)-1):
            if df_da.index[row_i] in cell_type_genes:
                if float(df_da.iloc[row_i,col_i])>0:
                    found_genes_cell.append(df_da.index[row_i])
                    print(f"match {df_da.index[row_i]} found")
        #print(found_genes_cell)
        if len(found_genes_cell)/len(cell_type_genes)>0.1:
            print("DA A9 found!")
            da_type_coord.append(int(col_i))
            #da_types_cell['DA_A9'].append(col_i)
    #da_types_cell[cell_type]=da_type_coord
    return(da_type_coord)


#da_types_cell={'DA_A9': [], 'DA_SNc' : [], 'DA_VTA1' : [], 'DA_VTA2' : [], 'DA_VTA3' : [], 'DA_VTA4' : []}

#da_types_cell = list(map(determine_cell_type, cell_types_da, df_da))
da_base = {k: determine_cell_type(v, df_da) for k, v in base_markers_da.items()}
da_types_cell = {k: determine_cell_type(v, df_da) for k, v in cell_types_da.items()}
l_max = max(map(len, da_types_cell.values()))
[row.extend([None]*(l_max - len(row))) for row in da_types_cell.values()]

da_types_cell
df_da_ann=pd.DataFrame.from_dict(da_types_cell)
df_da_ann

len(da_types_cell.values())
len(da_types_cell['DA_SNc'])
da_types_cell['DA_SNc']
len(da_types_cell['DA_A9'])

da_types_cell['DA_VTA4']
len(da_types_cell["DA_VTA1"])
len(da_types_cell['DA_VTA2'])
len(da_types_cell['DA_VTA3'])
len(da_types_cell['DA_VTA4'])
len(da_types_cell.items())
len(da_types_cell.keys())
#save the index information into a file. retrieve this information on R, take the column indeces
#and rename them according to the naming system

da_types_cell['DA_A9']

csv_p = '/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/'

df_da_ann.to_csv(csv_p + 'da_subtype_coord.csv', sep="\t")

#annotated_cells = list(map(iter_cell_types, cell_types_da))
da_types_cell
for type in cell_types_da.keys():
    cell_types_da[type]
#len(da_type_coord)
        #da_types_cell[type]=da_type_coord #save cell id into the given type

da_types_cell
da_types_cell["DA_A9"]
df_da


coord_list=[]

#Key: DA type, value: coords
da_type_coord = []

#get cells which contain markers for the given cell type. if cell has 60%+ of the markers
#of the given cell type, define it as this type. save info about which cell types correspond to the
#given marker

cell_types_da

da_type_coord = []
found_genes_cell=[]

da_types_cell={'DA_A9': [], 'DA_SNc' : [], 'DA_VTA1' : [], 'DA_VTA2' : [], 'DA_VTA4' : []}
for type in cell_types_da.keys():
    found_genes_cell=[]
    da_type_coord=[]
    for col_i, cell in enumerate(df_da_sub.columns):
        #start a gene counter - if 80% of the genes are expressed in the given group,
        #label them as the certain type
        print(f'Cell type atm: {cell_types_da[type]}')
        for row_i in range(len(df_da_sub.index)):
            if df_da_sub.index[row_i] in cell_types_da[type]:
                #print(f'da gene {df_da.index[row_i]} in data found...')
                if float(df_da_sub.iloc[row_i,col_i])>0:
                    found_genes_cell.append(df_da.index[row_i]) #the list saving the genes found per cell. we save the cell if X number of genes found
                    print(f"match {df_da_sub.index[row_i]} found")
                else:
                    continue
        #print(found_genes_cell)
        print(len(found_genes_cell)/len(cell_types_da[type]))
        if len(found_genes_cell)/len(cell_types_da[type])>0:
            #da_type_coord.append(col_i)
            #print(col_i)
            da_type_coord.append(col_i)
            #da_types_cell[type]=da_type_coord #save cell id into the given type
    da_types_cell[type].append(da_type_coord)
#############################################################

len(da_types_cell.values())
len(da_types_cell['DA_SNc'])
da_types_cell['DA_SNc']
da_types_cell['DA_A9']

da_types_cell['DA_VTA4']
len(da_types_cell["DA_VTA1"])
len(da_types_cell['DA_VTA2'])
len(da_types_cell['DA_VTA3'])
len(da_types_cell['DA_VTA4'])

da_types_cell['DA_SNc']


da_types_cell['DA_VTA4']
da_types_cell['DA_SNc']

ref.index = list(ref.index)
ref.index.iloc[4,0]


ref['cell type']
