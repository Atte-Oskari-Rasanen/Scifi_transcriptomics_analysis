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


from PIL import Image, ImageOps
import matplotlib.pyplot as plt
# Open the image form working directory


    

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
            split = img[i:i+tile_size, j:j+tile_size]
            tile_im=Image.fromarray(split)
            tile_im.save(tile_path + "/tile_" + str(j)+ "_" + str(i) + ".tif")

            img_tiles.append(split)
    return(img_tiles)

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


he_arr
cdna_arr
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

tile_2d_shape = (82, 82)
num_z = 1
num_r = 1
num_c = 2


synthetic_data = np.random.random(size=(num_r, num_c, num_z) + tile_2d_shape).astype(np.float32)
synthetic_data.shape


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
import os
tiles_arr[0]
tile_2d_shape = (82, 82)
num_z = 1
num_r = 1
num_c = 6#Write as a series of 3D tiffs.
import os, tempfile
from imageio import volread, volwrite

dir = tempfile.TemporaryDirectory()
dir

#saving the individual tiles
for r in range(num_r):
    for c in range(num_c):
        volwrite(os.path.join(dir.name, f"r{r}_c{c}.tiff"), tiles_arr[r, c])

import glob
print(glob.glob(str(dir).split(" ")[1].split(">")[0].split("'")[1].split("'")[0]+'/*'))

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

print(repr(stack))

#extracting the specific image from the numpy version (has the same dims as the stack ver) works fine but
#when transformed into stack, then taken xarray and then numpy via .values, then does not work. gives alwawys the 
#DAPI one only
stack.xarray[0,0,0]

stack.xarray[0,1,0].values

x_ar_DAPI=stack.xarray[0,1,0]
x_ar_DAPI.values
plt.imshow(x_ar_DAPI.values)
plt.show()
#stacks can be saved as the field of views

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
