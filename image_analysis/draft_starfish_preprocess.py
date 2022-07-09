from binascii import a2b_hex
from curses import window
from typing import Mapping, Union

import numpy as np
import starfish
from starfish.core.types import Coordinates, CoordinateValue, Axes

'''
we import the images as tiles and process them:
1. registration
2. autofluorescence removal (white top hat vs gaussian, maybe mexican hat e.g. to gene channel 2?)
3. finding the spots (blob log vs pixel wise)
4. pixel decode (hard vs soft)
Instead of finding spots to be decoded, it decodes every pixel and then connects
potential pixels with the same codeword into spots. 

the strength : works on dense data and noisy data where spot finding algorithms have a hard
time accurately detecting spots. 

the weakness : prone to false positives by decoding noise that would normally be ignored by
spot finding algorithms

'''
import os
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
# Open the image form working directory


'''
└── parent
    ├── slideA_1_1st_Cy3.5.TIF
    ├── slideA_1_1st_Cy3.TIF
    ├── slideA_1_1st_Cy5.TIF
    ├── slideA_1_1st_DAPI.TIF
    ├── slideA_1_1st_FITC.TIF
    ├── slideA_1_2nd_Cy3.5.TIF
    ├── slideA_1_2nd_Cy3.TIF
    ├── ...
------------------------------------
└── parent
    ├── <Fov1_name>
        └── <Fov1_name>
            ├── <target_name1>.tiff
            ├── ...
    ├── <Fov2_name>
        └── <Fov2_name>
            ├── <target_name1>.tiff
            ├── ...
-----------------------------------

Take each image tile, make into a single fov. So 1 Fov= 6 channels, 1, round, 1 z stack

'''    
base_path="/media/data/AtteR/projects/starfish/images/real_ims"

he_path=base_path+"/HE/MNM_509_D1_HE_Crop_coReg.png"

cdna_path=base_path+"/cDNA/MNM_509_D1_cDNA_Crop_coReg.ome.tif"

#genes_paths
c1_path=base_path+"/c1/c1_DARP32.tif"

c2_path=base_path+"/c2/c2_Chat.tif"

c3_path=base_path+"/c3/c3_PENK.tif"

c4_path=base_path+"/c4/c4_DRD1.tif"

im_paths=[he_path, cdna_path, c1_path, c2_path, c3_path, c4_path]


'''
before we put our data into spacetx format, we would process them into binaries right?

Our data structure: N_tiles FOVs, 1 rounds, 4 channels, and 1 z-planes
primary=4channels
dots=cdna
nuclei=he

N_tiles_x_fovs
prim_fovs = [
    [
        (0, 0, 0),
        (0, 1, 0),
        (0, 2, 0),
        (0, 3, 0),
        (0, 4, 0),
        (0, 5, 0),
    ]
]

he_fovs = [
    [
        (0, 0, 0),
    ]
]
'''


def start_points(size, split_size, overlap):
    points = [0]
    stride = int(split_size * (1-overlap))
    counter = 1
    while True:
        pt = stride * counter
        if pt + split_size >= size:
            points.append(size - split_size)
            break
        else:
            points.append(pt)
        counter += 1
    return points

import csv
def save_fov_coordinates(fov_path, fov_struc, coord_fovs_list):
    #save the coordinate file into each fov subdir
    with open(os.path.join(fov_path, "coordinates.csv"), "w") as fh:
        csv_writer = csv.DictWriter(
            fh,
            [
                'fov', 'round', 'ch', 'zplane',
                'xc_min', 'yc_min', 'zc_min', 'xc_max', 'yc_max', 'zc_max',
            ]
        )
        csv_writer.writeheader()
        for fov_id, (fov_info, coord_fov) in enumerate(zip(fov_struc, coord_fovs_list)):
            tile_coordinates = coord_fov.copy()
            tile_coordinates.update({
                'fov': fov_id,
                'round': fov_info[0],
                'ch': fov_info[1],
                'zplane': fov_info[2],
            })
            csv_writer.writerow(tile_coordinates)




def fovs_gen(im_arrs, tile_size, tile_points, base_path, fov_struc, image_type):
    #make this into a function that takes in the im_arrs list and iterates over the list
    #components
    #print(f'imported im_arrs shape: {im_arrs.shape}')
    
    print(f'imported im_arrs len: {len(im_arrs)}')

    coord_fovs={'xc_min':0, 'xc_max':0, 'yc_min':0, 'yc_max':0, 'zc_min':0.005, 'zc_max':0.010}
    '''
    determine the root dir (either nuclei or primary) inside which all the fovs are saved
    (equals to the number of tiles)
    '''
    if image_type=="nuclei":
            nucl_dir="nuclei"
            nucl_path = os.path.join(base_path, nucl_dir)
            try:
                os.mkdir(nucl_path)
            except FileExistsError:
                print("Directory already exists")
            base_path=nucl_path
    else:
        prim_dir="primary"
        prim_path = os.path.join(base_path, prim_dir)
        try:
            os.mkdir(prim_path)
        except FileExistsError:
            print("Directory already exists")
        base_path=prim_path
    #coord_fovs={'zc_min':0.005, 'zc_max':0.010}

    img_tiles=[]
    fov_n=0
    Y_points=tile_points[0]
    X_points=tile_points[1]
    for y_i, y in enumerate(tile_points[0][:-2]):
        for x_i, x in enumerate(tile_points[1][:-2]):
            coord_fovs={'zc_min':0.005, 'zc_max':0.010}

            '''
            Add here all the 6 channels, apply same procedure, name them, make their own
            fov subdir and save there. import images as lists. the HE saved into its own fov bits,
            the rest into the same fovs. make give coordinates as percs based on y and x points
            so:
            xmin=0
            xmax=xmin+tile_size
            '''

            #coordinates that will be saved into the coordinate file
            coord_fovs["yc_min"]=round(y/tile_points[0][-1], 6)
            coord_fovs["yc_max"]=round(tile_points[0][y_i+1]/tile_points[0][-1], 6)
            #coord_fovs["yc_max"]=Y_points[y+1]
            coord_fovs["xc_min"]=round(x/tile_points[1][-1], 6)
            #coord_fovs["xc_min"]=x
            coord_fovs["xc_max"]=round(X_points[x_i+1]/tile_points[1][-1], 6)
            coord_fovs_list=[]
            coord_fovs_list.append(coord_fovs)
            print(coord_fovs)
            fov_dir=f'fov_{fov_n}'
            fov_path = os.path.join(base_path, fov_dir)
            try:
                os.mkdir(fov_path)
                print("Directory '% s' created" % fov_path)
            except FileExistsError:
                print("Directory exists")

            #nuclei-f0-r2-c3-z33.tiff
            '''
            now you have the fov and fov coord dirs set up. in these places you save each tile
            dir structure
            └── Primary
                ├──<Fov_0>
                    ├── Primary-f0-r1-c1-z1.tiff
                    ├── Primary-f0-r1-c2-z1.tiff
                    ├── ...
                ├──<Fov_1>
                    ├── Primary-f1-r1-c1-z1.tiff
                    ├── Primary-f1-r1-c2-z1.tiff
                    ├── ...
            ------------------------------------
            └── nuclei
                ├──<Fov_0>
                    ├── nuclei-f0-r1-c1-z1.tiff
                    ├── nuclei-f0-r1-c2-z1.tiff
                    ├── ...
                ├──<Fov_1>
                    ├── nuclei-f1-r1-c1-z1.tiff
                    ├── nuclei-f1-r1-c2-z1.tiff
                    ├── ...
            '''
            
            #save the tile from each im to the fov dir
            if image_type=='nuclei':
                print(im_arrs[0].shape)
                print(type(im_arrs))

                split = im_arrs[0][y:y+tile_size, x:x+tile_size]
                tile_im=Image.fromarray(split)
                tile_im.save(f'{base_path}/{image_type}-f{fov_n}-r0-c0-z0.tif')
                img_tiles.append(split)
                xc_min_list.append(round(x/tile_points[1][-1], 6))
                yc_min_list.append(round(y/tile_points[0][-1]))
                zc_min_list.append(0.01)
                xc_max_list.append(round(X_points[x_i+1]/tile_points[1][-1],6))
                yc_max_list.append(round(Y_points[y_i+1]/tile_points[1][-1],6))
                zc_max_list.append(0.005)
                fov_list.append(fov_n)
                round_list.append(1)
                ch_list.append(1)
                zplane_list.append(1)
                
            else:
                for c,im in enumerate(im_arrs[1:]):
                    split = im[y:y+tile_size, x:x+tile_size]
                    tile_im=Image.fromarray(split)
                    tile_im.save(f'{base_path}/{image_type}-f{fov_n}-r0-c{c}-z0.tif')
                    img_tiles.append(split)
                    #fov,round,ch,zplane,xc_min,yc_min,zc_min,xc_max,yc_max,zc_max
                    xc_min_list.append(round(x/tile_points[1][-1], 6))
                    yc_min_list.append(round(y/tile_points[0][-1]))
                    zc_min_list.append(0.01)
                    xc_max_list.append(round(X_points[x_i+1]/tile_points[1][-1],6))
                    yc_max_list.append(round(Y_points[y_i+1]/tile_points[1][-1],6))
                    zc_max_list.append(0.005)
                    fov_list.append(fov_n)
                    round_list.append(1)
                    ch_list.append(c)
                    zplane_list.append(1)
            fov_n+=1
    df_coord = pd.DataFrame(list(zip(fov_list, round_list, ch_list, zplane_list, xc_min_list, yc_min_list,zc_min_list, xc_max_list, yc_max_list,zc_max_list)))
    df_coord.columns=['fov','round','ch','zplane','xc_min','yc_min','zc_min','xc_max','yc_max','zc_max']
    print("column names:")
    print(df_coord.columns)
    save_fov_coordinates(base_path, df_coord)
    return(img_tiles)

def tiles_gen(im_paths, tile_size,overlap, base_path):
    #first sort out the nuclei channel fovs, then the rest
    im_arrs=[]
    for path in im_paths:
        image = Image.open(path)
        #image = ImageOps.grayscale(image)
        img = np.asarray(image)
        x_pad=img.shape[0]%tile_size
        y_pad=img.shape[1]%tile_size
        img_padded = np.pad(img, ((x_pad,x_pad), (y_pad,y_pad)), constant_values=0.0, mode="constant")
        im_arrs.append(img_padded)
        print(f'Image imported, shape is {img_padded.shape}')
        #im_arrs[path.split("/")[8]]=img_padded
    print(f'im_arr length: {len(im_arrs)}')
    print(f'im_arr 1: {im_arrs[0].shape}')

    X_points = start_points(im_arrs[0].shape[0], tile_size, overlap)
    Y_points = start_points(im_arrs[0].shape[1], tile_size, overlap)
    tile_points=[Y_points, X_points]
    img_tiles=[]

    fov_structure_prim= [
            (0, 0, 0),
            (0, 1, 0),
            (0, 2, 0),
            (0, 3, 0),
            (0, 4, 0),
            (0, 5, 0),
        ]
    fov_structure_he= [
            (0, 0, 0)
        ]

    #this part was made in case i wanted to make a larger csv file containing all the possible fovs
    #N_fovs=len(im_arrs[0])/tile_size
    #prim_fovs = [[fov_structure_prim]*N_fovs]
    #nucl_fovs = [[fov_structure_he]*N_fovs]
    nuclei_arr=fovs_gen(im_arrs, tile_size, tile_points, base_path, fov_structure_he, "nuclei")
    #nuclei_arr=0
    prim_arr=fovs_gen(im_arrs[1:], tile_size, tile_points, base_path, fov_structure_prim, "primary")
    return([nuclei_arr, prim_arr])
    #pad images so that it is divisible by the tile size

base_path="/media/data/AtteR/projects/starfish/images/real_ims/FOVs"

all_arrs=tiles_gen(im_paths, 2000,0.1, base_path)


############################################################################
    

#only apply padding to the original one but iterate using the original shape as we fill in the
#values in the new version of the old array

#generate tiles
def gen_tiles(im_path,tile_path, tile_size, overlap):
    image = Image.open(im_path)
    #image = ImageOps.grayscale(image)
    img = np.asarray(image)

    img_h, img_w = img.shape

    def start_points(size, split_size, overlap):
        points = [0]
        stride = int(split_size * (1-overlap))
        counter = 1
        while True:
            pt = stride * counter
            if pt + split_size >= size:
                points.append(size - split_size)
                break
            else:
                points.append(pt)
            counter += 1
        return points


    X_points = start_points(img_w, tile_size, overlap)
    Y_points = start_points(img_h, tile_size, overlap)

    img_tiles=[]
    for i in Y_points:
        for j in X_points:
            #Add here all the 6 channels, apply same procedure, name them, make their own
            #fov subdir and save there. import images as lists
            split = img[i:i+tile_size, j:j+tile_size]
            tile_im=Image.fromarray(split)
            tile_im.save(tile_path + "/tile_" + str(j)+ "_" + str(i) + ".tif")
            img_tiles.append(split)
    return(img_tiles)

#make a function that takes in the X, Y coordinates

#class
# funct1: first import all images as np arrays, save into a list and return it
# funct2: takes one array from the list, doesnt matter which one as they all should be same size.
#it calculates the X and Y coordinates and iterates over them, making the split image of each array,
#creates a dict for the specific fov and saves the 6 different channels of the same tile into it


#make a class that has the function which takes in the x,y points and the list of 

def plot_image_grid(images, ncols=None, cmap='gray'):
    '''Plot a grid of images'''
    if not ncols:
        factors = [i for i in range(1, len(images)+1) if len(images) % i == 0]
        ncols = factors[len(factors) // 2] if len(factors) else len(images) // 4 + 1
    nrows = int(len(images) / ncols) + int(len(images) % ncols)
    imgs = [images[i] if len(images) > i else None for i in range(nrows * ncols)]
    f, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 2*nrows))
    axes = axes.flatten()[:len(imgs)]
    for img, ax in zip(imgs, axes.flatten()): 
        if np.any(img):
            if len(img.shape) > 2 and img.shape[2] == 1:
                img = img.squeeze()
            ax.imshow(img, cmap=cmap)
    


he_path="/media/data/AtteR/projects/starfish/images/HE.jpg"

#HE_arr=gen_tiles(he_path,82,0.5)
#plot_image_grid(HE_arr)
#plt.show()


Image.MAX_IMAGE_PIXELS = None


he_path="/media/data/AtteR/projects/starfish/images/real_ims/HE/MNM_509_D1_HE_Crop_coReg.png"
he_tiles_path="/media/data/AtteR/projects/starfish/images/real_ims/HE_tiles"

he_arr=gen_tiles(he_path, he_tiles_path,2000,0.1)


cdna_path="/media/data/AtteR/projects/starfish/images/real_ims/cDNA/MNM_509_D1_cDNA_Crop_coReg.ome.tif"
cdna_tiles_path="/media/data/AtteR/projects/starfish/images/real_ims/cDNA_tiles"
cdna_arr=gen_tiles(cdna_path, cdna_tiles_path,2000,0.1)


#genes_paths
c1_path="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes/unstacked/c1_DARP32.tif"
c1_path_tiles="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes_tiles/c1"
c1_arr=gen_tiles(c1_path, c1_path_tiles,2000,0.1)


c2_path="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes/unstacked/c2_Chat.tif"
c2_path_tiles="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes_tiles/c2"
c2_arr=gen_tiles(c2_path, c2_path_tiles,2000,0.1)

c3_path="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes/unstacked/c3_PENK.tif"
c3_path_tiles="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes_tiles/c3"
c3_arr=gen_tiles(c3_path, c3_path_tiles,2000,0.1)

c4_path="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes/unstacked/c4_DRD1.tif"
c4_path_tiles="/media/data/AtteR/projects/starfish/images/real_ims/RCA_genes_tiles/c4"
c4_arr=gen_tiles(c4_path, c4_path_tiles,2000,0.1)

'''
Maybe make the whole process more memory efficient by iteratively taking in an x number of tiles, then making
the starfish object and running analyses, then taking the rest and eventually merge the count matrices
'''


#make a script for segmenting them after they are part of starfish object or prior to it?


#list of the respective HE-DAPI pairs
resp_ims=[he_arr[100], cdna_arr[100], c1_arr[100], c2_arr[100], c3_arr[100], c4_arr[100]]
resp_ims[0].shape
#plt.imshow(resp_ims[0])
#plt.show()

tiles_arr = np.array([tile for tile in resp_ims])  #make list of the images into common array
tiles_arr=tiles_arr[:, np.newaxis, :, :]
tiles_arr=tiles_arr[np.newaxis,:]
tiles_arr.shape


'''
make tiles of example images, import as numpy arrays separately, combine?
'''
from typing import Mapping, Union

import numpy as np

from starfish.core.types import Coordinates, CoordinateValue, Axes



#Now build a FetchedTile and TileFetcher based on this data.
'''
The FetchedTile subclass defines the function you need for reading your images and the other
properties required by write_experiment_json() to construct slicedimage.Tiles

The TileFetcher subclass acts as the interface for write_experiment_json() to know 
where to get files to construct slicedimage.Tiles
'''
import functools
import skimage
print(skimage.__version__)

from starfish.experiment.builder import FetchedTile, TileFetcher
tile_2d_shape = (tiles_arr.shape[-2], tiles_arr.shape[-1])
tile_2d_shape
# We use this to cache images across tiles.  To avoid reopening and decoding the TIFF file, we use a
# single-element cache that maps between file_path and the array data.
@functools.lru_cache(maxsize=1)
def cached_read_fn(file_path) -> np.ndarray:
    return volread(file_path, format="tiff")

class DemoFetchedTile(FetchedTile):
    def __init__(self, filename, z, *args, **kwargs):
        self.filename = filename
        self.z = z

    @property
    def shape(self) -> Mapping[Axes, int]:
        return {
            Axes.Y: tile_2d_shape[0], Axes.X: tile_2d_shape[1],
        }

    @property
    def coordinates(self) -> Mapping[Union[str, Coordinates], CoordinateValue]:
        return {
            Coordinates.X: (0, 0.001),
            Coordinates.Y: (0, 0.001),
            Coordinates.Z: (0.001 * self.z, 0.001 * (self.z + 1)),
        }

    def tile_data(self) -> np.ndarray:
        return cached_read_fn(self.filename)[self.z]

class DemoTileFetcher(TileFetcher):
    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        return DemoFetchedTile(os.path.join(dir.name, f"r{r}_c{c}.tiff"), zplane_label)

#######################################

#Load the data as an ImageStack
'''
5-dimensional Image Tensor that labels each (z, y, x) tile with the round and channel and zplane it corresponds to.
ImageStacks can only be initialized with aligned Tilesets.
'''
from starfish import ImageStack
import os, tempfile
from imageio import volread, volwrite

num_z = 1
num_r = 1
num_c = 6#Write as a series of 3D tiffs.

#tile_2d_shape = (2000, 2000)
#synthetic_data = np.random.random(size=(num_r, num_c, num_z) + tile_2d_shape).astype(np.float32)
#synthetic_data.shape

dir = tempfile.TemporaryDirectory()
dir
#tiles_arr.shape

#saving the individual tiles
for r in range(num_r):
    for c in range(num_c):
        volwrite(os.path.join(dir.name, f"r{r}_c{c}.tiff"), tiles_arr[r, c])

import glob
print(glob.glob(str(dir).split(" ")[1].split(">")[0].split("'")[1].split("'")[0]+'/*'))

image = Image.open('/tmp/tmpr5bzdor_/r0_c4.tiff')
image = Image.open('/tmp/tmpr5bzdor_/r0_c0.tiff')
#image = ImageOps.grayscale(image)
img = np.asarray(image)
img.shape
plt.imshow(img)
plt.show()

#######################################

stack = ImageStack.from_tilefetcher(
    DemoTileFetcher(),
    {
        Axes.Y: tiles_arr.shape[3], Axes.X: tiles_arr.shape[4],
    },
    fov=0,
    rounds=range(tiles_arr.shape[0]),
    chs=range(tiles_arr.shape[1]),
    zplanes=range(tiles_arr.shape[2]),
    group_by=(Axes.ROUND, Axes.CH),
)
stack
####################################
'''
Example: 2 FOVs, 2 rounds, 1 channel, and 3 z-planes
'''
fovs = [
    [
        (0, 0, 0),
        (0, 0, 1),
        (0, 0, 2),
        (1, 0, 0),
        (1, 0, 1),
        (1, 0, 2),
    ],
    [
        (0, 0, 0),
        (0, 0, 1),
        (0, 0, 2),
        (1, 0, 0),
        (1, 0, 1),
        (1, 0, 2),
    ],
]



stack.xarray[0,0,0]
c3=stack.xarray[0,3,0].to_numpy
c3
tiles_arr.shape

stack.xarray.data
tiles_arr[0][0]
stack.xarray[0,1,0].values

x_ar_he=stack.xarray[0,1,0]
x_ar_he.values.shape
x_ar_cdna=stack.xarray[0,2,0]
x_ar_cdna.values.shape

plt.imshow(x_ar_cdna.values)
plt.show()
#stacks can be saved as the field of views

################
#fovs: prim. image and the auxiliary image
fovs[0]
for fov_id, fov in enumerate(fovs):
    for round_label, ch_label, zplane_label in fov:
        print(round_label)
        print(ch_label)
        print(zplane_label)
        print("======")

'''
primary images - genes
dots - cdna (contains all molecules of the experiment?)

rename the tiles:
<image_type>-f<fov_id>-r<round_label>-c<ch_label>-z<zplane_label>.

e.g:
nuclei-f0-r2-c3-z33.tiff

Because each image_type is treated as a separate set of images, you need a different
coordinates CSV file for each image_type. Therefore, each image_type must be converted
in its own directory containing all the 2D image files and and CSV file
'''
#extracting the specific image from the numpy version (has the same dims as the stack ver) works fine but
#when transformed into stack, then taken xarray and then numpy via .values, then does not work. gives alwawys the 
#DAPI one only


###################
from starfish import Codebook
sd = Codebook.synthetic_one_hot_codebook(n_round=1, n_channel=2, n_codes=2)
sd

'''
Build a 3-round 4-channel codebook where :code:`ACTA` is specified by intensity in round 0,
channel 1, and :code:`ACTB` is coded by fluorescence in channels 0, 1, and 2 of rounds 0,
1, and 2.
'''

import numpy as np
from starfish import Codebook

#make a numpy array of zeros, add pos. signals (1s) to parts where a signal for specific target should be found
data = np.zeros((1,2,1), dtype=np.uint8)
#data = np.zeros((1,2,1, 50,50), dtype=np.uint8)

data.shape
data

#so round 1, channel 2, z stack 1
data[0, 1, 0] = 1 # DAPI
data[0, 0, 0] = 1 # HE
from starfish.types import Axes, Features
from starfish import Codebook
codebook = [
    {
        Features.CODEWORD: [
            {Axes.ROUND.value: 0, Axes.CH.value: 1, Features.CODE_VALUE: 1},
            #{Axes.ROUND.value: 1, Axes.CH.value: 3, Features.CODE_VALUE: 1},
        ],
        Features.TARGET: "HE"
    },
    {
        Features.CODEWORD: [
            {Axes.ROUND.value: 0, Axes.CH.value: 2, Features.CODE_VALUE: 1},
            #{Axes.ROUND.value: 1, Axes.CH.value: 1, Features.CODE_VALUE: 1},
        ],
        Features.TARGET: "DAPI"
},
]
Codebook.from_code_array(codebook)

Codebook.from_numpy(['DAPI', 'HE'], n_channel=1, n_round=0, data=data)

###################

import matplotlib

from starfish.util.plot import imshow_plane
matplotlib.rcParams["figure.dpi"] = 150
f, (ax1, ax2, ax3) = plt.subplots(ncols=3)

# Plot first round and channel of projected ImageStack
imshow_plane(stack, sel={Axes.ROUND: 0, Axes.CH: 0, Axes.ZPLANE: 0}, ax=ax1, title='Z: 1')
imshow_plane(projection, sel={Axes.ROUND: 0, Axes.CH: 1}, ax=ax2, title='Max Projection')
# Plot ROI of projected image
selector = {Axes.CH: 0, Axes.ROUND: 0, Axes.X: (400, 600), Axes.Y: (550, 750)}
imshow_plane(projection, sel=selector, ax=ax3, title='Max Projection\nROI')



#register
from starfish.image import ApplyTransform, LearnTransform
from starfish.types import Axes

def register(imgs, dots, method = 'translation'):
    mip_imgs = imgs.reduce(dims = [Axes.CH, Axes.ZPLANE], func="max")
    mip_dots = dots.reduce(dims = [Axes.CH, Axes.ZPLANE], func="max")
    learn_translation = LearnTransform.Translation(reference_stack=mip_dots, axes=Axes.ROUND, upsampling=1000)
    transforms_list = learn_translation.run(mip_imgs)
    warp = ApplyTransform.Warp()
    registered_imgs = warp.run(imgs, transforms_list=transforms_list, in_place=False, verbose=True)
    return registered_imgs

#filter white top hat
from starfish.image import Filter
def filter_white_tophat(imgs, dots, masking_radius):
    wth = Filter.WhiteTophat(masking_radius=masking_radius)
    return wth.run(imgs), wth.run(dots)

# filter
imgs_wth, dots_wth = filter_white_tophat(imgs, dots, 15)

# register
registered_imgs = register(imgs_wth, dots_wth)

stack.get_image()

'''
All images can be accessed using get_image() with the name of the image type and any which 
rounds/chs/zplanes you want to include.
'''

from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json

outputdir = tempfile.TemporaryDirectory()

primary_tile_fetcher = PrimaryTileFetcher(primary_dir)
nuclei_tile_fetcher = NucleiTileFetcher(nuclei_dir)

# This is hardcoded for this example data set
primary_image_dimensions: Mapping[Union[str, Axes], int] = {
    Axes.ROUND: 2,
    Axes.CH: 1,
    Axes.ZPLANE: 3,
}
aux_images_dimensions: Mapping[str, Mapping[Union[str, Axes], int]] = {
    "nuclei": {
        Axes.ROUND: 2,
        Axes.CH: 1,
        Axes.ZPLANE: 3,
    },
}

write_experiment_json(
    path=outputdir.name,
    fov_count=2,
    tile_format=ImageFormat.TIFF,
    primary_image_dimensions=primary_image_dimensions,
    aux_name_to_dimensions=aux_images_dimensions,
    primary_tile_fetcher=primary_tile_fetcher,
    aux_tile_fetcher={"nuclei": nuclei_tile_fetcher},
    dimension_order=(Axes.ROUND, Axes.CH, Axes.ZPLANE)
)
