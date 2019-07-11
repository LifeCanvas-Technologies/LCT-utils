"""
A module to enable command-line cropping of TIFF files and stacks
"""

import argparse
import glob
import sys
import os
import platform
from typing import List, Tuple, Sequence, Union
import numpy
import tifffile
from tqdm import tqdm


def parse_args(args=sys.argv[1:]):
    """Reads in command line arguments to script"""
    parser = argparse.ArgumentParser(description="Crop tiff files.")
    parser.add_argument("--input",
                        required=True,
                        help="A tiff file to crop.")
    parser.add_argument("--output",
                        help="The path to the output tiff file (if not "
                             "specified, will return dimensions of input tiff).")
    parser.add_argument("--crop-x",
                        help="Specify as \"--crop-x [x-min],[x-max]\".")
    parser.add_argument("--crop-y",
                        help="Specify as \"--crop-y [y-min],[y-max]\".")
    parser.add_argument("--crop-z",
                        help="Specify as \"--crop-z [z-min],[z-max]\".")
    parser.add_argument("--split-output",
                        help="Specify whether or not the output should be saved as a collection of tiffs"
                             "Note: the OUTPUT should be a path to a folder, not a file.",
                        action="store_true")

    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=sys.argv[1:]):
    args = parse_args(args)
    if args.output is None:
        check_dims(input_glob_or_list=args.input,
                   print_out=True)
    else:
        crop_tiff(input_glob_or_list=args.input,
                  output_path=args.output,
                  crop_x=args.crop_x,
                  crop_y=args.crop_y,
                  crop_z=args.crop_z,
                  split_output=args.split_output)


def crop_tiff(input_glob_or_list: Union[str, List[str]],
              output_path: str,
              crop_x: Union[str, List[int], Tuple[int, int]] = None,
              crop_y: Union[str, List[int], Tuple[int, int]] = None,
              crop_z: Union[str, List[int], Tuple[int, int]] = None,
              split_output: bool = False) -> None:
    """Main function for cropping tiffs

    Parameters
    ----------
    input_glob_or_list : Union[str, List[str]]
        String containing input glob or list of file paths to individual tiffs
    output_path : str
        Output file path, shoudl either be a directory or tiff file, depending on the value of split_output
    crop_x : Union[str, List[int], Tuple[int, int]]
        String in format "[x_min,x_max]" or Sequence of ints in format [x_min, x_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff
    crop_y : Union[str, List[int], Tuple[int, int]]
        String in format "[y_min,y_max]" or Sequence of ints in format [y_min, y_max] specifying the minimum
        and maximum y-values (inclusive) to include in the outputted cropped tiff
    crop_z : Union[str, List[int], Tuple[int, int]]
        String in format "[z_min,z_max]" or Sequence of ints in format [z_min, z_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff
    split_output : bool
        True if output should be a directory of 2D tiffs, False if output should be single multi-page tiff file

    Returns
    -------
    None
    """

    tiff_paths = glob_to_list(input_glob_or_list)

    # Ensures TIFF stack is outputted into a directory
    if split_output:
        if platform.system() == "Windows":
            if output_path[-1] != "\\":
                output_path += "\\"
        elif platform.system() == "Linux" or platform.system == "Darwin":
            if output_path[-1] != "/":
                output_path += "/"

    x_dim, y_dim, z_dim = check_dims(input_glob_or_list)
    crop_x, crop_y, crop_z = crop_str_to_tuple(x_dim, y_dim, z_dim, crop_x, crop_y, crop_z)

    if len(tiff_paths) == 1:
        tiff_ndarray_reshape = reshape_single_tiff(tiff_paths[0])
        crop_single_tiff(tiff_ndarray_reshape,
                         output_filepath=output_path,
                         output_tiff_stack=split_output,
                         crop_x=crop_x,
                         crop_y=crop_y,
                         crop_z=crop_z)
    else:
        crop_multi_tiff(tiff_paths,
                        output_filepath=output_path,
                        output_tiff_stack=split_output,
                        crop_x=crop_x,
                        crop_y=crop_y,
                        crop_z=crop_z)


def check_dims(input_glob_or_list, print_out: bool = True) -> Tuple[int, int, int]:
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


def crop_str_to_tuple(x_dim: int = None,
                      y_dim: int = None,
                      z_dim: int = None,
                      crop_x: Union[str, List[int], Tuple[int, int]] = None,
                      crop_y: Union[str, List[int], Tuple[int, int]] = None,
                      crop_z: Union[str, List[int], Tuple[int, int]] = None) \
        -> Tuple[List[int], List[int], List[int]]:
    """Turns two comma separated ints in string form to a list of two int
    Expects dimensions of the original tiff file to use as defaults in the event that cropping dimensions
    are not specified. If crop_x, crop_y, and crop_z are guaranteed to be specified, x_dim, y_dim, and z_dim
    are not used.

    Parameters
    ----------
    x_dim : int
        Number of voxels in the x-direction in the original tiff
    y_dim : int
        Number of voxels in the y-direction in the original tiff
    z_dim : int
        Number of voxels in the z-direction in the original tiff
    crop_x : Union[str, List[int], Tuple[int, int]]
        String in format "[x_min,x_max]" or Sequence of ints in format [x_min, x_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff
    crop_y : Union[str, List[int], Tuple[int, int]]
        String in format "[y_min,y_max]" or Sequence of ints in format [y_min, y_max] specifying the minimum
        and maximum y-values (inclusive) to include in the outputted cropped tiff
    crop_z : Union[str, List[int], Tuple[int, int]]
        String in format "[z_min,z_max]" or Sequence of ints in format [z_min, z_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff

    Returns
    -------
    tuple[List[int], List[int], List[int]]
        Tuple of 3 lists, each of 2 ints, in format ( [x_min, x_max], [y_min, y_max], [z_min, z_max] )
    """


    if x_dim is None and crop_x is None:
        raise ValueError("Need to specify at least x_dim or crop_x")
    if y_dim is None and crop_y is None:
        raise ValueError("Need to specify at least y_dim or crop_y")
    if z_dim is None and crop_z is None:
        raise ValueError("Need to specify at least z_dim or crop_z")

    if x_dim is None:
        x_dim = float("Inf")
    if y_dim is None:
        y_dim = float("Inf")
    if z_dim is None:
        z_dim = float("Inf")

    if crop_x is None:
        crop_x_list = [1, x_dim]
    elif type(crop_x) is list:
        crop_x_list = crop_x
    elif type(crop_x) is tuple:
        crop_x_list = list(crop_x)
    else:
        crop_x_list = [int(s) for s in crop_x.split(",")]
    if crop_x_list[0] < 0:
        raise ValueError("Minimum x-cropping value must be positive")
    if crop_x_list[0] > crop_x_list[1]:
        raise ValueError("Minimum x-cropping value must be smaller"
                         "than maximum x-cropping value")
    if crop_x_list[1] > x_dim:
        raise ValueError("Maximum x-cropping value must be smaller than"
                         "the dimensions of the image")
    if len(crop_x_list) != 2:
        raise ValueError("2 components need to be specified, not {}!".format(len(crop_x)))

    if crop_y is None:
        crop_y_list = [1, y_dim]
    elif type(crop_y) is list:
        crop_y_list = crop_y
    elif type(crop_y) is tuple:
        crop_y_list = list(crop_y)
    else:
        crop_y_list = [int(s) for s in crop_y.split(",")]
    if crop_y_list[0] < 0:
        raise ValueError("Minimum y-cropping value must be positive")
    if crop_y_list[0] > crop_y_list[1]:
        raise ValueError("Minimum y-cropping value must be smaller"
                         "than maximum y-cropping value")
    if crop_y_list[1] > y_dim:
        raise ValueError("Maximum y-cropping value must be smaller than"
                         "the dimensions of the image")
    if len(crop_y_list) != 2:
        raise ValueError("2 components need to be specified, not {}!".format(len(crop_y)))

    if crop_z is None:
        crop_z_list = [1, z_dim]
    elif type(crop_z) is list:
        crop_z_list = crop_z
    elif type(crop_z) is tuple:
        crop_z_list = list(crop_z)
    else:
        crop_z_list = [int(s) for s in crop_z.split(",")]
    if crop_z_list[0] < 0:
        raise ValueError("Minimum z-cropping value must be positive")
    if crop_z_list[0] > crop_z_list[1]:
        raise ValueError("Minimum z-cropping value must be smaller"
                         "than maximum z-cropping value")
    if crop_z_list[1] > z_dim:
        raise ValueError("Maximum z-cropping value must be smaller than"
                         "the dimensions of the image")
    if len(crop_z_list) != 2:
        raise ValueError("2 components need to be specified, not {}!".format(len(crop_z)))

    return crop_x_list, crop_y_list, crop_z_list


def reshape_single_tiff(file_path: str) -> numpy.ndarray:
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
        tiff_ndarray_reshape = tiff_ndarray[numpy.newaxis, ...]
    elif tiff_ndarray.ndim == 3:
        if tiff_ndarray.shape[-1] < 10:
            # Interpret as 2D TIFF with extra channels
            tiff_ndarray_reshape = tiff_ndarray[numpy.newaxis, ...]
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


def get_single_tiff_dims(tiff_ndarray: numpy.ndarray):
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


def crop_single_tiff(
        tiff_ndarray: numpy.ndarray,
        output_filepath: str,
        output_tiff_stack: bool,
        crop_x: Sequence[int],
        crop_y: Sequence[int],
        crop_z: Sequence[int]) -> None:
    """Crops a TIFF ndarray, saving as a TIFF file

    Parameters
    ----------
    tiff_ndarray : numpy.ndarray
        A numpy ndarray, already padded as needed to include 3 or 4
        dimensions, depending on whether or not there are extra channels,
        in the following order: [z, x, y, <extra channels>]
    output_filepath : str
        A file path or directory to which the resulting TIFF file will be saved. If the
        file or directory does not exist, it will be created.
    output_tiff_stack : bool
        True if the cropped tiff should be saved as a collection of 2D tiffs,
        False if it should be saved in a multipage tiff
    crop_x : Sequence[int]
        A list or tuple in the format [x_min, x_max] such that the resulting TIFF
        file is cropped to include all x-dimension pixels in the range
        (x_min, x_max), inclusive, of the original TIFF file
    crop_y : Sequence[int]
        A list or tuple in the format [y_min, y_max] such that the resulting TIFF
        file is cropped to include all y-dimension pixels in the range
        (y_min, y_max), inclusive, of the original TIFF file
    crop_z : Sequence[int]
        A list or tuple in the format [z_min, z_max] such that the resulting TIFF
        file is cropped to include all z-dimension pixels in the range
        (z_min, z_max), inclusive, of the original TIFF file

    Returns
    -------
    None

    """
    x_min, x_max = crop_x
    y_min, y_max = crop_y
    z_min, z_max = crop_z

    output_tiff_size = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
    if output_tiff_size > 2 * 1024 ** 3:
        bigtiff = True
    else:
        bigtiff = False

    tiff_ndarray_cropped = tiff_ndarray[z_min:z_max + 1, y_min:y_max + 1, x_min:x_max + 1]

    if output_tiff_stack is False:
        with tifffile.TiffWriter(output_filepath, bigtiff=bigtiff) as file:
            file.save(tiff_ndarray_cropped)
    else:
        if not os.path.exists(output_filepath):
            print("Creating directory at {}.".format(output_filepath))
            os.mkdir(output_filepath)
        tiff_num = 0
        for z in tqdm(range(z_min - 1, z_max)):
            tiff_num_str = str(tiff_num)
            while len(tiff_num_str) < numpy.log10(z_max - z_min + 1) + 1:
                tiff_num_str = "0" + tiff_num_str
            tiff_filename = output_filepath + "img_" + tiff_num_str + ".tiff"
            with tifffile.TiffWriter(tiff_filename, bigtiff=bigtiff) as file:
                file.save(tiff_ndarray_cropped[z])
            tiff_num += 1


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


def crop_multi_tiff(
        tiff_paths: List[str],
        output_filepath: str,
        output_tiff_stack: bool = False,
        crop_x: Sequence[int] = None,
        crop_y: Sequence[int] = None,
        crop_z: Sequence[int] = None) -> None:
    """Crops TIFF stacks without loading the entire stack into memory at once
    Note: assumes at least one TIFF in the stack fits in memory

    Parameters
    ----------
    tiff_paths : List[str]
        A list of file path to the tiffs in the stack, in order of ascending z-coordinate
    output_filepath : str
        A file path or directory to which the resulting TIFF file will be saved. If the
        file or directory does not exist, it will be created.
    output_tiff_stack : bool
        True if the cropped tiff should be saved as a collection of 2D tiffs,
        False if it should be saved in a multipage tiff
    crop_x : Sequence[int]
        A list or tuple in the format [x_min, x_max] such that the resulting TIFF
        file is cropped to include all x-dimension pixels in the range
        (x_min, x_max), inclusive, of the original TIFF file
    crop_y : Sequence[int]
        A list or tuple in the format [y_min, y_max] such that the resulting TIFF
        file is cropped to include all y-dimension pixels in the range
        (y_min, y_max), inclusive, of the original TIFF file
    crop_z : Sequence[int]
        A list or tuple in the format [z_min, z_max] such that the resulting TIFF
        file is cropped to include all z-dimension pixels in the range
        (z_min, z_max), inclusive, of the original TIFF file

    Returns
    -------
    None
    """
    x_min, x_max = crop_x
    y_min, y_max = crop_y
    z_min, z_max = crop_z

    output_tiff_size = (x_max - x_min)*(y_max - y_min)*(z_max - z_min)
    if output_tiff_size > 2 * 1024**3:
        bigtiff = True
    else:
        bigtiff = False

    if output_tiff_stack is False:
        with tifffile.TiffWriter(output_filepath, bigtiff=bigtiff) as file:
            for z in tqdm(range(z_min - 1, z_max)):
                tiff_ndarray = tifffile.imread(tiff_paths[z])
                tiff_ndarray_cropped = tiff_ndarray[y_min:y_max, x_min:x_max]
                file.save(tiff_ndarray_cropped)
    else:
        if not os.path.exists(output_filepath):
            print("Creating directory at {}.".format(output_filepath))
            os.mkdir(output_filepath)
        tiff_num = 0
        for z in tqdm(range(z_min - 1, z_max)):
            tiff_num_str = str(tiff_num)
            while len(tiff_num_str) < numpy.log10(z_max - z_min + 1) + 1:
                tiff_num_str = "0" + tiff_num_str
            tiff_filename = output_filepath + "img_" + tiff_num_str + ".tiff"
            with tifffile.TiffWriter(tiff_filename, bigtiff=bigtiff) as file:
                tiff_ndarray = tifffile.imread(tiff_paths[z])
                tiff_ndarray_cropped = tiff_ndarray[y_min:y_max + 1, x_min:x_max + 1]
                file.save(tiff_ndarray_cropped)
            tiff_num += 1


if __name__ == "__main__":
    main()
