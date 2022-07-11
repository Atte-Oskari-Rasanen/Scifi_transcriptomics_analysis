import numpy as np
import os
from PIL import Image, ImageOps
Image.MAX_IMAGE_PIXELS = None
import json

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

def fovs_gen(im_arrs, tile_size, tile_points, base_path, manifest_file, image_type):
   
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
            FOV_dir= {"default_tile_format": "TIFF", "dimensions": ["z","xc","x","yc","y","zc","c","r"],"extras": {},
            "shape": {
                "c": 1,
                "r": 0,
                "z": 0
            },
            "tiles": [],
            "version": "0.1.0"
            }


    else:
        prim_dir="primary"
        prim_path = os.path.join(base_path, prim_dir)
        try:
            os.mkdir(prim_path)
        except FileExistsError:
            print("Directory already exists")
        base_path=prim_path
        FOV_dir= {"default_tile_format": "TIFF", "dimensions": ["z","xc","x","yc","y","zc","c","r"],"extras": {},
            "shape": {
                "c": 5,
                "r": 0,
                "z": 0
            },
            "tiles": [],
            "version": "0.1.0"
            }

    #coord_fovs={'zc_min':0.005, 'zc_max':0.010}
    FOV_dir_orig=FOV_dir.copy()
    Y_points=tile_points[0]
    X_points=tile_points[1]
    fov_i=0
    for y_i, y in enumerate(tile_points[0][:-2]):
        for x_i, x in enumerate(tile_points[1][:-2]):

            #coordinate_tiles=[]
            #save the tile from each im to the fov dir
            if image_type=='nuclei':
                split = im_arrs[0][y:y+tile_size, x:x+tile_size]                
                tile_im=Image.fromarray(split)
                tilename=f'{image_type}-f{fov_i}-r0-c0-z0.tiff'
                tile_im.save(f'{base_path}/{tilename}')

                xmin=round(x/tile_points[1][-1], 6)
                xmax=round(X_points[x_i+1]/im_arrs[0].shape[0],6)
                ymin=round(y/tile_points[0][-1])
                ymax=round(Y_points[y_i+1]/im_arrs[0].shape[1],6)
                zmin=0.01
                zmax=0.005
                tilename=f'{image_type}-f{fov_i}-r0-c0-z0.tiff'
                fov_i_dir={'coordinates': {'xc': [xmin, xmax], 'yc': [ymin, ymax],
                'zc': [zmin, zmax]}, 'file': tilename,
                'indices': {'c': 0, 'r': 0, 'z': 0}, 'tile_format': 'TIFF',
                'tile_shape': {'x': tile_size, 'y': tile_size}}
                FOV_dir["tiles"].append(fov_i_dir)
                fov_key=f'fov_{fov_i}'
                im_name=f'{base_path}/{image_type}-images-fov_{fov_i}.json'
                current_fov_dir={fov_key: im_name}
                manifest_file["contents"].update(current_fov_dir)

                with open(f'{base_path}/{image_type}-images-fov_{fov_i}.json', "w") as jsonfile:
                    json.dump(FOV_dir,jsonfile)

            else:
                for c,im in enumerate(im_arrs[1:]):
                    split = im[y:y+tile_size, x:x+tile_size]                
                    tile_im=Image.fromarray(split)
                    tilename=f'{image_type}-f{fov_i}-r0-c{c}-z0.tiff'

                    tile_im.save(f'{base_path}/{tilename}')
                    #img_tiles.append(split)
                    #fov,round,ch,zplane,xc_min,yc_min,zc_min,xc_max,yc_max,zc_max
                    xmin=round(x/tile_points[1][-1], 6)
                    xmax=round(X_points[x_i+1]/im_arrs[0].shape[0],6)
                    ymin=round(y/tile_points[0][-1])
                    ymax=round(Y_points[y_i+1]/im_arrs[0].shape[1],6)
                    zmin=0.01
                    zmax=0.005
                    fov_i_dir={'coordinates': {'xc': [xmin, xmax], 'yc': [ymin, ymax],
                    'zc': [zmin, zmax]}, 'file': tilename,
                    'indices': {'c': c, 'r': 0, 'z': 0}, 'tile_format': 'TIFF',
                    'tile_shape': {'x': tile_size, 'y': tile_size}}
                    FOV_dir["tiles"].append(fov_i_dir)
                with open(f'{base_path}/{image_type}-images-fov_{fov_i}.json', "w") as jsonfile:
                    json.dump(FOV_dir,jsonfile)

                FOV_dir["tiles"]=[]
                fov_key=f'fov_{fov_i}'
                im_name=f'{base_path}/{image_type}-images-fov_{fov_i}.json'
                current_fov_dir={fov_key: im_name}
                manifest_file["contents"].update(current_fov_dir)
            FOV_dir=FOV_dir_orig
            fov_i+=1

    with open(f'{base_path}/{image_type}_images.json', 'w') as jsonfile:
        json.dump(manifest_file,jsonfile) 
    print(f'manifest file saved with all the fovs as: {base_path}/{image_type}_images.json')
    return(FOV_dir)


def tiles_gen(im_paths, tile_size,overlap, base_path,exp_file, codebook, manifest_file):
    #first sort out the nuclei channel fovs, then the rest
    im_arrs=[]
    for path in im_paths:
        image = Image.open(path)
        #image = ImageOps.grayscale(image)
        img = np.asarray(image)
        x_pad=img.shape[0]%tile_size
        y_pad=img.shape[1]%tile_size
        #img_padded = np.pad(img, ((x_pad,x_pad), (y_pad,y_pad)), constant_values=0.0, mode="constant")
        im_arrs.append(img)
        print(f'Image imported, shape is {img.shape}')

    X_points = start_points(im_arrs[0].shape[0], tile_size, overlap)
    Y_points = start_points(im_arrs[0].shape[1], tile_size, overlap)
    tile_points=[Y_points, X_points]

    nuclei_arr=fovs_gen(im_arrs, tile_size, tile_points, base_path,manifest_file,"nuclei")
    #nuclei_arr=0
    prim_arr=fovs_gen(im_arrs[1:], tile_size, tile_points, base_path,manifest_file,"primary")
    #prim_arr=0

    #pass the manifest files into the functions and update as you go through the loops
    with open(f'{base_path}/experiment.json', "w") as jsonfile:
        json.dump(exp_file,jsonfile) 
    with open(f'{base_path}/codebook.json', "w") as jsonfile:
        json.dump(codebook,jsonfile) 

    print(f'experiment file saved with all the fovs as: {base_path}/experiment.json')

    return(im_arrs)
    #pad images so that it is divisible by the tile size


