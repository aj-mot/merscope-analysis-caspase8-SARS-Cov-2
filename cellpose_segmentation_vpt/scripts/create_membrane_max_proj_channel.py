#Combines membrane images from MERSCOPE data directory and creates a maximum projection image for each z level
#combines this with dapi and saves it as a tiff file for each z level
#Image will be saved as mosaic_Cellbound_MAX_z{z}.tif, where z is the z level

import os
import argparse
from aicsimageio import AICSImage
from tifffile import TiffWriter
from joblib import Parallel, delayed
import dask.array as da

def save_tif(img_combined, save_path, idx):
    img_combined = img_combined.compute()
    with TiffWriter(os.path.join(save_path, f'mosaic_Cellbound_MAX_z{idx}.tif'), bigtiff=True) as tif:
        options = dict(
            tile=(512, 512),
            compression='zlib',
            photometric='minisblack',
            metadata={
                'axes': 'YX',
                'PhysicalSizeX': 0.108,
                'PhysicalSizeY': 0.108
            }
        )
        tif.write(img_combined, **options)

def process_images(img_path, save_path):
    assert os.path.exists(img_path), f"Image path {img_path} does not exist."
    os.makedirs(save_path, exist_ok=True)

    dapi_list = []
    membrane_1 = []
    membrane_2 = []
    membrane_3 = []

    for file in os.listdir(img_path):
        if "dapi" in file.lower():
            dapi_list.append(os.path.join(img_path, file))
        elif "cellbound1" in file.lower():
            membrane_1.append(os.path.join(img_path, file))
        elif "cellbound2" in file.lower():
            membrane_2.append(os.path.join(img_path, file))
        elif "cellbound3" in file.lower():
            membrane_3.append(os.path.join(img_path, file))
    #check if any of membrane lists are empty
    assert len(membrane_1) == len(membrane_2) == len(membrane_3), "Membrane images have different lengths."
    for idx in range(len(membrane_3)):
        print(f"Processing z level {idx}")
        membrane1_img = da.squeeze(AICSImage(membrane_1[idx]).dask_data)
        membrane2_img = da.squeeze(AICSImage(membrane_2[idx]).dask_data)
        membrane3_img = da.squeeze(AICSImage(membrane_3[idx]).dask_data)

        channel_max = da.maximum(membrane1_img, membrane2_img, membrane3_img)
        save_tif(channel_max, save_path, idx)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process membrane images and save combined TIFF files.")
    parser.add_argument("--img_paths", nargs='+', type=str, help="List of paths to the input images.")
    parser.add_argument("--save_paths", nargs='+', type=str, help="List of paths to save the output images.")
    args = parser.parse_args()

    # Ensure both img_paths and save_paths are provided and have the same length
    assert args.img_paths is not None, "Please provide the paths to the input images."
    assert args.save_paths is not None, "Please provide the paths to save the output images."
    assert len(args.img_paths) == len(args.save_paths), "The number of input paths must match the number of output paths."

    # Function to process a each pair of img_path and save_path
    def process_pair(img_path, save_path):
        assert os.path.exists(img_path), f"Image path {img_path} does not exist."
        process_images(img_path, save_path)
        print(f"Finished processing images from {img_path} and saved to {save_path}.")

    # Use joblib to parallelize the processing
    Parallel(n_jobs=-1)(delayed(process_pair)(img_path, save_path) for img_path, save_path in zip(args.img_paths, args.save_paths))
