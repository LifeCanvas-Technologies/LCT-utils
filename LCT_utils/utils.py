"""Utility functions for LCT-utils
"""

import glob
import sys
from typing import Tuple, List, Union

import numpy as np
import tifffile
import tqdm.auto as tqdm


def glob_to_list(glob_or_list: Union[List[str], Tuple[str, ...], str]) -> List[str]:
    """Converts globs to sorted lists of file paths
    If a list or tuple is passed in, it will be sorted and returned as a list

    Parameters
    ----------
    glob_or_list : Union[List[str], Tuple[str, ...], str]
        If a glob is passed in, a corresponding list of sorted file paths will be returned, else if a list or
        tuple is passed in, it will be sorted and returned as a list

    Returns
    -------
    List[str]
        A sorted list of file paths

    """
    if type(glob_or_list) is list:
        output_list = sorted(glob_or_list)
        if len(output_list) < 1:
            raise Exception("Input list of length 0!")
    elif type(glob_or_list) is tuple:
        output_list = sorted(list(glob_or_list))
        if len(output_list) < 1:
            raise Exception("Input tuple of length 0!")
    elif type(glob_or_list) is str:
        output_list = sorted(glob.glob(glob_or_list))
        if len(output_list) < 1:
            raise Exception("No files found at {}".format(glob_or_list))
    else:
        raise TypeError("Input was not a glob, list, or tuple")

    return output_list


def get_tiff_dims(input_glob_or_list: Union[List[str], Tuple[str, ...], str],
                  print_out: bool = True) -> Tuple[int, int, int]:
    """Checks the dimensions of the tiff file(s) at the given glob

    Parameters
    ----------
    input_glob_or_list : Union[str, List[str], Tuple[str]]
        Either a UNIX-style glob corresponding to the respective tiff files or a list/tuple of file paths
    print_out : bool
        True if the dimensions of the tiff should be printed to the console, False otherwise

    Returns
    -------
    Tuple[int, int, int]
        Tuple of ints corresponding to dimensions in the x, y, and z directions respectively

    """
    tiff_paths = glob_to_list(input_glob_or_list)

    if len(tiff_paths) == 1:
        tiff_ndarray = reshape_single_tiff(tiff_paths[0])
        x_dim, y_dim, z_dim = get_single_tiff_dims(tiff_ndarray)
    else:
        x_dim, y_dim, z_dim = get_multi_tiff_dims(tiff_paths)

    if print_out:
        dims_str = str(x_dim) + ", " + str(y_dim) + ", " + str(z_dim)
        print("Dimensions (x, y, z) (width, height, depth) [units: voxels]  : ", dims_str)

    return x_dim, y_dim, z_dim


def reshape_single_tiff(file_path: str) -> np.ndarray:
    """Given a file path, returns reformatted ndarray and its dimensions

    Note: assumes there are no more than 9 extra channels and that
    there are at least 10 pixels in the y dimension

    Parameters
    ----------
    file_path : str
        A path to a tiff file, e.g. "/path-to-tiffs/example.tiff"

    Returns
    -------
    numpy.ndarray
        A reformatted ndarray with dimensions
        in the following order: [z, y, x, <extra channels>], where extra
        channels are optional (including RGB, CMYK, alpha, etc.)
    """
    tiff_ndarray = tifffile.imread(file_path)

    if tiff_ndarray.ndim == 2:
        # Interpret as 2D greyscale TIFF
        tiff_ndarray_reshape = tiff_ndarray[np.newaxis, ...]
    elif tiff_ndarray.ndim == 3:
        if tiff_ndarray.shape[-1] < 10:
            # Interpret as 2D TIFF with extra channels
            tiff_ndarray_reshape = tiff_ndarray[np.newaxis, ...]
        else:
            # Interpret as 3D TIFF
            tiff_ndarray_reshape = tiff_ndarray.copy()
    elif tiff_ndarray.ndim == 4:
        # Interpret as 3D TIFF with extra channels
        tiff_ndarray_reshape = tiff_ndarray.copy()
    else:
        sys.stderr.write("Not yet configured to handle this kind of TIFF.")
        sys.exit(1)

    return tiff_ndarray_reshape


def get_single_tiff_dims(tiff_ndarray: np.ndarray):
    """Returns dimensions of a tiff represented by a reshaped array

    Parameters
    ----------
    tiff_ndarray : numpy.ndarray
        A reshaped numpy array represented the tiff at hand

    Returns
    -------
    x_dim : int
        The number of voxels in the x-direction
    y_dim : int
        The number of voxels in the y-direction
    z_dim : int
        The number of voxels in the z-direction
    """
    z_dim, y_dim, x_dim = tiff_ndarray.shape[0:3]
    return x_dim, y_dim, z_dim


def get_multi_tiff_dims(tiff_paths: List[str]) -> Tuple[int, int, int]:
    """Given a UNIX-style glob, returns dimensions of the TIFF stack and a list of file paths

    Parameters
    ----------
    tiff_paths : List[str]
        a list of file paths to each individual tiff in the tiff stack

    Returns
    -------
    x_dim : int
        The number of voxels in the x-direction
    y_dim : int
        The number of voxels in the y-direction
    z_dim : int
        The number of voxels in the z-direction
    """
    z_dim = len(tiff_paths)
    y_dim, x_dim = tifffile.imread(tiff_paths[0]).shape[0:2]

    return x_dim, y_dim, z_dim


def parse_crop_dims_str(crop_x: Union[str, List[float], Tuple[float, float]] = None,
                        crop_y: Union[str, List[float], Tuple[float, float]] = None,
                        crop_z: Union[str, List[float], Tuple[float, float]] = None) \
        -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """Turns two comma separated ints or floats in string form to a tuple of two ints or floats
    In the event that no cropping is specified in a given dimension, will return bounds of -Inf and +Inf

    Parameters
    ----------
    crop_x : Union[str, List[float], Tuple[float, float]]
        String in format "[x_min,x_max]" or Sequence of ints in format [x_min, x_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff
    crop_y : Union[str, List[float], Tuple[float, float]]
        String in format "[y_min,y_max]" or Sequence of ints in format [y_min, y_max] specifying the minimum
        and maximum y-values (inclusive) to include in the outputted cropped tiff
    crop_z : Union[str, List[float], Tuple[float, float]]
        String in format "[z_min,z_max]" or Sequence of ints in format [z_min, z_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff

    Returns
    -------
    Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]
        Tuple of 3 tuples, each of 2 ints or floats (depending on inputs),
        in format ( (x_min, x_max), (y_min, y_max), (z_min, z_max) )
    """
    # X-axis cropping
    if crop_x is None:
        crop_x_tuple = (-float("inf"), float("inf"))
    elif type(crop_x) is list:
        crop_x_tuple = tuple(crop_x)
    elif type(crop_x) is tuple:
        crop_x_tuple = crop_x
    else:
        try:
            crop_x_tuple = tuple([int(s) for s in crop_x.split(",")])
        except ValueError:
            crop_x_tuple = tuple([float(s) for s in crop_x.split(",")])
    if crop_x_tuple[0] > crop_x_tuple[1]:
        raise ValueError("Minimum x-cropping value must be smaller"
                         "than maximum x-cropping value")
    if len(crop_x_tuple) != 2:
        raise ValueError("2 components need to be specified, not {}!".format(len(crop_x)))

    # Y-axis cropping
    if crop_y is None:
        crop_y_tuple = (-float("inf"), float("inf"))
    elif type(crop_y) is list:
        crop_y_tuple = tuple(crop_y)
    elif type(crop_y) is tuple:
        crop_y_tuple = crop_y
    else:
        try:
            crop_y_tuple = tuple([int(s) for s in crop_y.split(",")])
        except ValueError:
            crop_y_tuple = tuple([float(s) for s in crop_y.split(",")])
    if crop_y_tuple[0] > crop_y_tuple[1]:
        raise ValueError("Minimum y-cropping value must be smaller"
                         "than maximum y-cropping value")
    if len(crop_y_tuple) != 2:
        raise ValueError("2 components need to be specified, not {}!".format(len(crop_y)))

    # Z-axis cropping
    if crop_z is None:
        crop_z_tuple = (-float("inf"), float("inf"))
    elif type(crop_z) is list:
        crop_z_tuple = tuple(crop_z)
    elif type(crop_z) is tuple:
        crop_z_tuple = crop_z
    else:
        try:
            crop_z_tuple = tuple([int(s) for s in crop_z.split(",")])
        except ValueError:
            crop_z_tuple = tuple([float(s) for s in crop_z.split(",")])
    if crop_z_tuple[0] > crop_z_tuple[1]:
        raise ValueError("Minimum z-cropping value must be smaller"
                         "than maximum z-cropping value")
    if len(crop_z_tuple) != 2:
        raise ValueError("2 components need to be specified, not {}!".format(len(crop_z)))

    return crop_x_tuple, crop_y_tuple, crop_z_tuple


def validate_tiff_crop_dims(x_dim: int,
                            y_dim: int,
                            z_dim: int,
                            crop_x: Tuple[float, float],
                            crop_y: Tuple[float, float],
                            crop_z: Tuple[float, float]) \
        -> Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]:
    """Ensures cropping dimensions are valid for a tiff, attempts to fix if not

    Parameters
    ----------
    x_dim : int
        The number of voxels in the x-dimension
    y_dim : int
        The number of voxels in the y-direction
    z_dim : int
        The number of voxels in the z-direction
    crop_x : Tuple[float, float]
        A tuple pair specifying the minimum and maximum bounds to include in the x-direction after the crop
    crop_y : Tuple[float, float]
        A tuple pair specifying the minimum and maximum bounds to include in the y-direction after the crop
    crop_z : Tuple[float, float]
        A tuple pair specifying the minimum and maximum bounds to include in the z-direction after the crop

    Returns
    -------
    Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]
        A tuple triple of validated tuple pairs in the following format:
        ((x_min, x_max), (y_min, y_max), (z_min, z_max))
    """
    # Handle out of range bounds
    # NOTE: must happen first to handle Inf and -Inf properly
    if crop_x[0] < 1:
        # warn("x-cropping lower bound less than 1, setting to 1")
        crop_x = (1, crop_x[1])
    if crop_x[1] > x_dim:
        # warn("x-cropping upper bound exceeds TIFF bounds, setting to {}".format(x_dim))
        crop_x = (crop_x[0], x_dim)
    if crop_y[0] < 1:
        # warn("y-cropping lower bound less than 1, setting to 1")
        crop_y = (1, crop_y[1])
    if crop_y[1] > y_dim:
        # warn("y-cropping upper bound exceeds TIFF bounds, setting to {}".format(y_dim))
        crop_y = (crop_y[0], y_dim)
    if crop_z[0] < 1:
        # warn("z-cropping lower bound less than 1, setting to 1")
        crop_z = (1, crop_z[1])
    if crop_z[1] > z_dim:
        # warn("z-cropping upper bound exceeds TIFF bounds, setting to {}".format(z_dim))
        crop_z = (crop_z[0], z_dim)

    # Handle float bounds
    if type(crop_x[0]) is float:
        # warn("x-cropping lower bound is float, converting to int")
        crop_x = (int(crop_x[0]), crop_x[1])
    if type(crop_x[1]) is float:
        # warn("x-cropping upper bound is float, converting to int")
        crop_x = (crop_x[0], int(crop_x[1]))
    if type(crop_y[0]) is float:
        # warn("y-cropping lower bound is float, converting to int")
        crop_y = (int(crop_y[0]), crop_y[1])
    if type(crop_y[1]) is float:
        # warn("y-cropping upper bound is float, converting to int")
        crop_y = (crop_y[0], int(crop_y[1]))
    if type(crop_z[0]) is float:
        # warn("z-cropping lower bound is float, converting to int")
        crop_z = (int(crop_z[0]), crop_z[1])
    if type(crop_z[1]) is float:
        # warn("z-cropping upper bound is float, converting to int")
        crop_z = (crop_z[0], int(crop_z[1]))

    return crop_x, crop_y, crop_z


def pretty_print_json_coords(coords: List[List[float]], output_path: str) -> None:
    """A quick function to improve formatting and speed of writing out JSON coordinate triples
    NOTE: can only handle coordinate triples at this time.

    Parameters
    ----------
    coords : List[List[float]]
        List of lists, where each inner list specifies a coordinate triple
    output_path : str
        Path to which the resulting json will be outputted

    Returns
    -------
    None

    """
    # Ensure coords formatting is okay
    try:
        __, __, __ = coords[0]
    except ValueError:
        print("Cannot yet handle coords in this format!")
        raise

    with open(output_path, "w") as fd:
        first_time = True
        fd.write("[\n")
        for x, y, z in tqdm.tqdm(coords):
            if first_time:
                first_time = False
            else:
                fd.write(",\n")
            fd.write("  [%.1f,%.1f,%.1f]" % (x, y, z))
        fd.write("]\n")
