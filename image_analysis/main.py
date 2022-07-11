from binascii import a2b_hex
from tkinter import image_types
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
os.chdir('/media/data/AtteR/projects/starfish')
from images_reformat_spacetx import *
# Open the image form working directory


'''
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
The function fovs_gen takes in a dict which follows the following structure:
----------------
FOV_dir structure
----------------
Main keys:
- default_tile_format - x
- dimensions - x
- extras - x
- shape - x 
- tiles - [fov_tile_coordinate_1, fov_tile_coordinate_2]
- version

-----------------------------
fov_tile_coordinate structure:
-----------------------------
{'coordinates': {'xc': [0.0, 0.112], 'yc': [0.0, 0.0539],
'zc': [0.0, 0.0001]}, 'file': 'primary-fov_000-c1-r0-z0.tiff',
'indices': {'c': 1, 'r': 0, 'z': 0}, 'tile_format': 'TIFF',
'tile_shape': {'x': 1120, 'y': 539}

'''

manifest_file={"version": "0.0.0",
    "contents": {},
    "extras": 'none'
    }
exp_file={
"version": "5.0.0",
"images": {
    "primary": "primary/primary_images.json",
    "nuclei": "nuclei/nuclei_images.json"

},
"codebook": "codebook.json",
"extras": {
    "is_space_tx_cool": 'True'
    }
}    

#0-HE, 1-cDNA, 2-DARP32, 3-Chat, 4-PENK, 5-DRD1

codebook={
  "version": "0.0.0",
  "mappings": [
    {
      "codeword": [
        {"c": 0, "r": 0, "v": 1},
        {"c": 2, "r": 0, "v": 1}

      ],
      "target": "DARP32"
    },
    {
      "codeword": [
        {"c": 0, "r": 0, "v": 1},
        {"c": 3, "r": 0, "v": 1}

      ],
      "target": "Chat"
    },
    {
      "codeword": [
        {"c": 0, "r": 0, "v": 1},
        {"c": 4, "r": 0, "v": 1}

      ],
      "target":"PENK"
    },

    {
      "codeword": [
        {"c": 0, "r": 0, "v": 1},
        {"c": 5, "r": 0, "v": 1}

      ],
      "target":"DRD1"
    }
  ]
}


base_path="/media/data/AtteR/projects/starfish/images/real_ims/FOVs"

all_arrs=tiles_gen(im_paths, 2000,0.1, base_path, exp_file, codebook, manifest_file)


##################################################################
from starfish import Experiment
#AttributeError: 'list' object has no attribute 'decode' --- maybe issues with how the fov files were created? 
#you appended to the tiles 

'''
Starfish loads data by referencing the top-level experiment.json objects in a SpaceTx Format
dataset. The main way to load data on your machine is through the Experiment constructor
as follows:
'''

#make a function for visualising all 6 channel tiles next to each other and then another one 
#for visualising an individual one

import xarray as xr

#could put the HE as the nuclei image... and maybe cdna as the dots?
experiment = Experiment.from_json(base_path +"/experiment.json")
experiment

'''
fov_0: <starfish.FieldOfView>
  Primary Image: <slicedimage.TileSet (c: 4, z: 1, r: 3, x: 2000, y: 2000)>
'''
e_fov1=experiment['fov_0'].get_image('primary')
e_fov1=experiment['fov_0'].get_image('nuclei')



e_fov2=experiment['fov_2'].get_image('primary')
a=e_fov2.xarray
np_arr=np.squeeze(a.to_numpy())

def show_im(e_fov, c):
    a=e_fov.xarray[0,c,0]
    np_arr=np.squeeze(a.to_numpy())
    plt.imshow(np_arr)
    plt.show()

def compare_ims(e_fov):
    a=e_fov.xarray
    np_arr=np.squeeze(a.to_numpy())
    f, axs = plt.subplots(2, 3, figsize=(15, 15))
    axs=axs.ravel()
    for c in range(np_arr.shape[0]):
        axs[c].imshow(np_arr[c])
    plt.show()

compare_ims(e_fov2)
show_im(e_fov2,4)

from starfish.types import Axes

projected_imgs = e_fov1.reduce({Axes.CH}, func="max")
print(projected_imgs)
show_im(projected_imgs,0)

plt.imshow(projected_imgs)
import matplotlib
from starfish.util.plot import diagnose_registration

matplotlib.rcParams["figure.dpi"] = 250
diagnose_registration(projected_imgs, {Axes.CH:0}, {Axes.CH:1}, {Axes.CH:2}, {Axes.CH:3})

experiment.codebook
#get the total

############################################################################
from starfish.image import ApplyTransform, LearnTransform
from starfish.types import Axes

#The maximum intensity projection of any round from primary images can also be used in lieu of a dots image

def register(imgs, dots, method = 'translation'):
    mip_imgs = imgs.reduce(dims = [Axes.CH, Axes.ZPLANE], func="max")
    mip_dots = dots.reduce(dims = [Axes.CH, Axes.ZPLANE], func="max")
    learn_translation = LearnTransform.Translation(reference_stack=mip_dots, axes=Axes.ROUND, upsampling=1000)
    transforms_list = learn_translation.run(mip_imgs)
    warp = ApplyTransform.Warp()
    registered_imgs = warp.run(imgs, transforms_list=transforms_list, in_place=False, verbose=True)
    return registered_imgs
register(e_fov1,dots)
from starfish import Codebook
sd = Codebook.synthetic_one_hot_codebook(n_round=1, n_channel=6, n_codes=4)
sd

'''
Build a 3-round 4-channel codebook where :code:`ACTA` is specified by intensity in round 0,
channel 1, and :code:`ACTB` is coded by fluorescence in channels 0, 1, and 2 of rounds 0,
1, and 2.
'''


import numpy as np

'''
#make a numpy array of zeros, add pos. signals (1s) to parts where a signal for specific target should be found
data = np.zeros((1,2,1), dtype=np.uint8)
#data = np.zeros((1,2,1, 50,50), dtype=np.uint8)

data.shape
data


The size of the codeword arrays should equal the “r” and “c” dimensions of the primary imagestack


#so round 1, channel 2, z stack 1
data[0, 1, 0] = 1 # HE
data[0, 2, 0] = 1 # cDNA
data[0, 3, 0] = 1 #DARP32
data[0, 4, 0] = 1 #Chat
data[0, 5, 0] = 1 #PENK
data[0, 6, 0] = 1 #DRD1
'''
from starfish.types import Axes, Features
from starfish import Codebook
codebook = [
    {
        Features.CODEWORD: [
            {Axes.ROUND.value: 0, Axes.CH.value: 1, Features.CODE_VALUE: 1},
            {Axes.ROUND.value: 0, Axes.CH.value: 3, Features.CODE_VALUE: 1},
        ],
        Features.TARGET: "DARP32"
    },
    {
        Features.CODEWORD: [
            {Axes.ROUND.value: 0, Axes.CH.value: 1, Features.CODE_VALUE: 1},
            {Axes.ROUND.value: 1, Axes.CH.value: 4, Features.CODE_VALUE: 1},
        ],
        Features.TARGET: "Chat"
},
    {
        Features.CODEWORD: [
            {Axes.ROUND.value: 0, Axes.CH.value: 1, Features.CODE_VALUE: 1},
            {Axes.ROUND.value: 1, Axes.CH.value: 5, Features.CODE_VALUE: 1},
        ],
        Features.TARGET: "PENK"
},
    {
        Features.CODEWORD: [
            {Axes.ROUND.value: 0, Axes.CH.value: 1, Features.CODE_VALUE: 1},
            {Axes.ROUND.value: 1, Axes.CH.value: 6, Features.CODE_VALUE: 1},
        ],
        Features.TARGET: "DRD1"
},

]
Codebook.from_code_array(codebook)

Codebook.from_numpy(['DARP32', 'Chat', 'PENK', 'DRD1'], n_channel=6, n_round=0, data=data)


)


############################################################################
import functools
from imageio import volread
from skimage.io import imread
from typing import Mapping, Union

from starfish.experiment.builder import FetchedTile
from starfish.types import Axes, Coordinates


# a 2D read function
def read_fn(file_path) -> np.ndarray:
    return imread(file_path)


# example of a cached 3D read function
# not used in this example
@functools.lru_cache(maxsize=1)
def cached_3D_read_fn(file_path) -> np.ndarray:
    return volread(file_path)


# subclass FetchedTile
class RNATile(FetchedTile):

    def __init__(
            self,
            file_path: str,
            coordinates: Mapping[Union[str, Coordinates], tuple]
    ) -> None:
        """Parser for a tile.

        Parameters
        ----------
        file_path : str
            location of the tiff
        coordinates : Mapping[Union[str, Coordinates], tuple]
            the coordinates for the selected tile, extracted from the metadata
        """
        self.file_path = file_path

        # coordinates must match shape
        self._coordinates = coordinates

    @property
    def shape(self) -> Mapping[Axes, int]:
        return {Axes.Y: 10, Axes.X: 10}  # hard coded for this example

    @property
    def coordinates(self):
        return self._coordinates

    def tile_data(self) -> np.ndarray:
        return read_fn(self.file_path)

#############
# physical coordinates for two FOVs

coordinates_of_fovs = [
    {
        Coordinates.X: (0.0, 0.1),
        Coordinates.Y: (0.0, 0.1),
        Coordinates.Z: (0.005, 0.010),
    },
    {
        Coordinates.X: (0.1, 0.2),
        Coordinates.Y: (0.0, 0.1),
        Coordinates.Z: (0.005, 0.010),
    },
]
coordinates_of_fovs=pd.read_csv('/media/data/AtteR/projects/starfish/images/real_ims/FOVs/nuclei/coordinates.csv',sep='\t')
import pandas as pd
coord=pd.read_csv('/media/data/AtteR/projects/starfish/images/real_ims/FOVs/primary/coordinates.csv',sep='\t')
coordinates_of_fovs = coord.loc[:, ~coord.columns.str.contains('^Unnamed')]

prim_path="/media/data/AtteR/projects/starfish/images/real_ims/FOVs/primary"
for subdir, dirs, files in os.walk(prim_path):
    for file in files:
        if file=='coordinates.csv':
            print(os.path.join(subdir, file))


from starfish.experiment.builder import TileFetcher
'''
how can the coordinate file be same for the primary and nuclei images? they consist of different
number of channels. and the codebook cant be the same for both either as codebook structure is 
N(rounds)xN(channels)
'''
class PrimaryTileFetcher(TileFetcher):

    def __init__(self, input_dir: str) -> None:
        self.input_dir = os.path.join(input_dir)
        self.num_z = 1

    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = f"primary-f{fov_id}-r{round_label}-c{ch_label}-z{zplane_label}.tiff"
        return RNATile(os.path.join(self.input_dir, filename), coordinates_of_fovs[fov_id])




#give each fov dir shape: c:6, r:1, z:1


FOV_dir["tiles"][0]
FOV_dir["tiles"][1]['file']

#so for each coordinate, save th

class NucleiTileFetcher(TileFetcher):

    def __init__(self, input_dir: str) -> None:
        self.input_dir = os.path.join(input_dir)
        self.num_z = 1

    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = f"nuclei-f{fov_id}-r{round_label}-c{ch_label}-z{zplane_label}.tiff"
        return RNATile(os.path.join(self.input_dir, filename), coordinates_of_fovs[fov_id])
#############
from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json

import tempfile

outputdir = tempfile.TemporaryDirectory()

primary_tile_fetcher = PrimaryTileFetcher(primary_dir)

nuclei_tile_fetcher = NucleiTileFetcher("/media/data/AtteR/projects/starfish/images/real_ims/FOVs/nuclei")

# This is hardcoded for this example data set
primary_image_dimensions: Mapping[Union[str, Axes], int] = {
    Axes.ROUND: 1,
    Axes.CH: 6,
    Axes.ZPLANE: 1,
}
aux_images_dimensions: Mapping[str, Mapping[Union[str, Axes], int]] = {
    "nuclei": {
        Axes.ROUND: 1,
        Axes.CH: 1,
        Axes.ZPLANE: 1,
    },
}
'''
data must be in spacetx format if you want to load your experiment into a starfish pipeline.
spacetx uses json files to organise single.plane tiffs.

'''

#whats our field of views? this varies based on the number of channels?

write_experiment_json(
    path=outputdir.name,
    fov_count=169,    
    tile_format=ImageFormat.TIFF,
    primary_image_dimensions=primary_image_dimensions,
    aux_name_to_dimensions=aux_images_dimensions,
    primary_tile_fetcher=primary_tile_fetcher,
    aux_tile_fetcher={"nuclei": nuclei_tile_fetcher},
    dimension_order=(Axes.ROUND, Axes.CH, Axes.ZPLANE)
)


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
