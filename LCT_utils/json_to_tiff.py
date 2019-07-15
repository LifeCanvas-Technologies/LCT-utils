"""
A module to transform json file of coordinates to TIFF files (assumes coordinates are in x, y, z format)
"""

import sys
import os
import json
import argparse
import multiprocessing as mp
from typing import List

import numpy as np
import tifffile
import tqdm.auto as tqdm

try:
    from .utils import get_tiff_dims
except ModuleNotFoundError:
    from LCT_utils.utils import get_tiff_dims


def parse_args(args=sys.argv[1:]):
    """Reads in command line arguments to script"""
    parser = argparse.ArgumentParser(description="Turn JSON coordinates into a TIFF stack")

    parser.add_argument("--input",
                        required=True,
                        help="A tiff file to crop.")
    parser.add_argument("--output",
                        required=True,
                        help="The path to the output directory of tiffs")
    parser.add_argument("--reference-tiff",
                        required=True,
                        help="A UNIX-style glob to the TIFF stack from which the coordinates are derived")

    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=sys.argv[1:]):
    args = parse_args(args)
    x_dim, y_dim, z_dim = get_tiff_dims(args.reference_tiff)

    with open(args.input) as fd:
        full_coords = json.load(fd)

    coords_by_z = [[] for __ in range(z_dim)]

    # Rounds z to nearest integer
    for x, y, z in full_coords:
        coords_by_z[np.intc(z) - 1].append([np.intc(x), np.intc(y)])

    # Create output directory if it does not yet exist
    if not os.path.exists(args.output):
        print("Creating directory at {}.".format(args.output))
        os.mkdir(args.output)

    with mp.Pool() as pool:
        results = []
        for z in range(z_dim):
            filename = args.output + "img_" + str(z).zfill(int(np.floor(np.log10(z_dim)) + 1)) + ".tiff"
            results.append(pool.apply_async(write_tiff, (z, coords_by_z, filename, x_dim, y_dim)))

        [result.get() for result in tqdm.tqdm(results)]


def write_tiff(z: int,
               coords_by_z: List[List[int, int, int]],
               filename: str,
               x_dim,
               y_dim):
    """Write a plane in the output TIFF stack

    Parameters
    ----------
    z : int
        z of current plane
    coords_by_z : List
        List of Lists, each inner list is coordinate triple
    filename : str
        File path to TIFF to be written
    x_dim : int
        x-dimension of a plane
    y_dim : int
        y-dimension of a plane

    Returns
    -------
    None
    """
    dtype = np.uint8
    dtype_max = np.iinfo(dtype).max

    plane = np.ones((y_dim, x_dim), dtype=dtype)
    x_coords = [xy[0] for xy in coords_by_z[z]]
    y_coords = [xy[1] for xy in coords_by_z[z]]
    plane[y_coords, x_coords] = dtype_max - 1

    with tifffile.TiffWriter(filename) as file:
        file.save(plane)


if __name__ == "__main__":
    main()