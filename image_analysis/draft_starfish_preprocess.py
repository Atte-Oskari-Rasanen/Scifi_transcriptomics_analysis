from binascii import a2b_hex
from curses import window
from typing import Mapping, Union

import numpy as np

from starfish.core.types import Coordinates, CoordinateValue, Axes

import cv2


from PIL import Image, ImageOps
import matplotlib.pyplot as plt
# Open the image form working directory
def make_arr(im_path, tile_size):

    image = Image.open(im_path)
    image = ImageOps.grayscale(image)
    data = np.asarray(image)
    im = Image.fromarray(data)

    def split(img, window_size, margin):
        sh = list(img.shape)
        sh[0], sh[1] = sh[0] + margin * 2, sh[1] + margin * 2
        img_ = np.zeros(shape=sh)
        img_[margin:-margin, margin:-margin] = img
        stride = window_size
        step = window_size + 2 * margin
        nrows, ncols = img.shape[0] // window_size, img.shape[1] // window_size
        splitted = []
        for i in range(nrows):
            for j in range(ncols):
                h_start = j*stride
                v_start = i*stride
                cropped = img_[v_start:v_start+step, h_start:h_start+step]
                splitted.append(cropped)
        tiles_arr = np.array([tile for tile in splitted])  #make list of the images into common array
        return tiles_arr

    #make the script take into account that if the window_size larger than either x or y of the array, then give notification!
    out = split(data, window_size=tile_size, margin=1)
    out.shape


    out[0].shape
    image2 = Image.fromarray(out[0])
    #plt.imshow(image2, interpolation='nearest')

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
        
    plot_image_grid(out)
    plt.show()
    return(out)
he_path="/media/data/AtteR/projects/starfish/images/HE.jpg"
dapi_path="/media/data/AtteR/projects/starfish/images/DAPI.jpg"
HE_arr=make_arr(he_path, 80)
DAPI_arr=make_arr(dapi_path,80)

HE_arr[0].shape
DAPI_arr[0].shape

#list of the respective HE-DAPI pairs
resp_ims=[HE_arr[0], DAPI_arr[0]]
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

#Write as a series of 3D tiffs.
import os, tempfile
from imageio import volread, volwrite

dir = tempfile.TemporaryDirectory()

for r in range(num_r):
    for c in range(num_c):
        volwrite(os.path.join(dir.name, f"r{r}_c{c}.tiff"), tiles_arr[r, c])

import glob
print(glob.glob(str(dir).split(" ")[1].split(">")[0].split("'")[1].split("'")[0]+'/*'))

#######################################

#Now build a FetchedTile and TileFetcher based on this data.
'''
The FetchedTile subclass defines the function you need for reading your images and the other
properties required by write_experiment_json() to construct slicedimage.Tiles

The TileFetcher subclass acts as the interface for write_experiment_json() to know 
where to get files to construct slicedimage.Tiles
'''
import functools
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

stack = ImageStack.from_tilefetcher(
    DemoTileFetcher(),
    {
        Axes.Y: tile_2d_shape[0], Axes.X: tile_2d_shape[1],
    },
    fov=0,
    rounds=range(num_r),
    chs=range(num_c),
    zplanes=range(num_z),
    group_by=(Axes.ROUND, Axes.CH),
)
print(repr(stack))

from starfish import Codebook
sd = Codebook.synthetic_one_hot_codebook(n_round=1, n_channel=2, n_codes=2)
sd


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
