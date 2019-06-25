"""
A module to enable command-line cropping of TIFF files and stacks
"""

import argparse
import glob
import sys
import os
import platform
from typing import List
import numpy
import tifffile
from tqdm import tqdm


def _parse_args(args=sys.argv[1:]):
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
                            "Note: the OUTPUT should be tje path to a folder, not a file.",
                        action="store_true")

    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=sys.argv[1:]):
    args = _parse_args(args)

    # Ensures TIFF stack is outputted into a directory
    output_filepath = args.output
    if args.split_output:
        if platform.system() == "Windows":
            if output_filepath[-1] != "\\":
                output_filepath += "\\"
        elif platform.system() == "Linux" or platform.system == "Darwin":
            if output_filepath != "/":
                output_filepath += "/"

    if "*" in args.input:
        x_dim, y_dim, z_dim, tiff_paths = get_tiff_stack_dims(args.input)
        if args.output is None:
            if z_dim == 1:
                units_str = "[units: pixels]"
            else:
                units_str = "[units: voxels]"
            dims_str = str(x_dim) + ", " + str(y_dim) + ", " + str(z_dim)
            print("Dimensions (x, y, z) (width, height, depth)", units_str, " : ", dims_str)
        else:
            crop_tiff_stack(tiff_paths,
                            x_dim,
                            y_dim,
                            z_dim,
                            output_filepath=args.output,
                            output_tiff_stack=args.split_output,
                            crop_x=args.crop_x,
                            crop_y=args.crop_y,
                            crop_z=args.crop_z)
    else:
        x_dim, y_dim, z_dim, tiff_ndarray_reshape = reshape_and_get_tiff_dims(args.input)

        if args.output is None:
            if z_dim == 1:
                units_str = "[units: pixels]"
            else:
                units_str = "[units: voxels]"
            dims_str = str(x_dim) + ", " + str(y_dim) + ", " + str(z_dim)
            print("Dimensions (x, y, z) (width, height, depth)", units_str, " : ", dims_str)
        else:
            crop_tiffs(tiff_ndarray_reshape,
                       x_dim=x_dim,
                       y_dim=y_dim,
                       z_dim=z_dim,
                       output_filepath=args.output,
                       output_tiff_stack=args.split_output,
                       crop_x=args.crop_x,
                       crop_y=args.crop_y,
                       crop_z=args.crop_z)


def reshape_and_get_tiff_dims(file_path: str) -> numpy.ndarray:
    """Given a file path, returns reformatted ndarray and its dimensions

    Note: assumes there are no more than 9 extra channels and that
    there are at least 10 pixels in the y dimension

    Parameters
    ----------
    file_path : str
        A path to a tiff file, e.g. "/path-to-tiffs/example.tiff"

    Returns
    -------
    {
        "x_dim" : int,
        "y_dim" : int,
        "z_dim" : int,
        "tiff_ndarray_reshape": numpy.ndarray
    }
        "x_dim", "y_dim", "z_dim" are ints with the number of pixels in
        the x, y, and z (width, height, depth) dimensions, respectively

        "tiff_ndarray_reshape" is a reformatted ndarray with dimensions
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
        if tiff_ndarray.shape[-1] < 10:
            # Interpret as 3D TIFF with extra channels
            tiff_ndarray_reshape = tiff_ndarray.copy()
        else:
            sys.stderr.write("Not yet configured to handle this kind of TIFF.")
            sys.exit(1)
    else:
        sys.stderr.write("Not yet configured to handle this kind of TIFF.")
        sys.exit(1)

    z_dim, y_dim, x_dim = tiff_ndarray_reshape.shape[0:3]
    return x_dim, y_dim, z_dim, tiff_ndarray_reshape


def crop_tiffs(
        tiff_ndarray: numpy.ndarray,
        x_dim: int,
        y_dim: int,
        z_dim: int,
        output_filepath: str,
        output_tiff_stack: bool,
        crop_x: str = None,
        crop_y: str = None,
        crop_z: str = None) -> bool:
    """Crops a TIFF ndarray, saving as a TIFF file
    Parameters
    ----------
    tiff_ndarray : numpy.ndarray
        A numpy ndarray, already padded as needed to include 3 or 4
        dimensions, depending on whether or not there are extra channels,
        in the following order: [z, x, y, <extra channels>]
    x_dim : int
        The number of pixels or voxels in the x-dimension
    y_dim : int
        The number of pixels or voxels in the y-dimension
    z_dim : int
        The number of pixels or voxels in the z-dimension
    output_filepath : str
        A file path or directory to which the resulting TIFF file will be saved. If the
        file or directory does not exist, it will be created.
    output_tiff_stack : bool
        True if the cropped tiff should be saved as a collection of 2D tiffs,
        False if it should be saved in a multipage tiff
    crop_x : str
        A string in the format "[x_min],[x_max]" such that the resulting TIFF
        file is cropped to include all x-dimension pixels in the range
        (x_min, x_max), inclusive, of the original TIFF file
    crop_y : str
        A string in the format "[y_min],[y_max]" such that the resulting TIFF
        file is cropped to include all y-dimension pixels in the range
        (y_min, y_max), inclusive, of the original TIFF file
    crop_z : str
        A string in the format "[z_min],[z_max]" such that the resulting TIFF
        file is cropped to include all z-dimension pixels in the range
        (z_min, z_max), inclusive, of the original TIFF file

    Returns
    -------
    bool
        True if successful, False otherwise

    """
    if crop_x is None:
        x_min = 1
        x_max = x_dim
    else:
        x_min, x_max = [int(s) for s in crop_x.split(",")]
        if x_min < 0:
            raise ValueError("Minimum x-cropping value must be positive")
        if x_min > x_max:
            raise ValueError("Minimum x-cropping value must be smaller",
                             "than maximum x-cropping value")
        if x_max > x_dim:
            raise ValueError("Maximum x-cropping value must be smaller than",
                             "the dimensions of the image")

    if crop_y is None:
        y_min = 1
        y_max = y_dim
    else:
        y_min, y_max = [int(s) for s in crop_y.split(",")]
        if y_min < 0:
            raise ValueError("Minimum y-cropping value must be positive")
        if y_min > y_max:
            raise ValueError("Minimum y-cropping value must be smaller",
                             "than maxmimum y-cropping value")
        if y_max > y_dim:
            raise ValueError("Maximum y-cropping value must be smaller than",
                             "the dimensions of the image")

    if crop_z is None:
        z_min = 1
        z_max = z_dim
    else:
        z_min, z_max = [int(s) for s in crop_z.split(",")]
        if z_min < 0:
            raise ValueError("Minimum z-cropping value must be positive")
        if z_min > z_max:
            raise ValueError("Minimum z-cropping value must be smaller",
                             "than maxmimum z-cropping value")
        if z_max > z_dim:
            raise ValueError("Maximum z-cropping value must be smaller than",
                             "the dimensions of the image")

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

    return True


def get_tiff_stack_dims(file_path_glob: str):
    """Given a UNIX-style glob, returns dimensions of the TIFF stack and a list of file paths

    :param file_path_glob:
    :return:
    """
    tiff_paths = sorted(glob.glob(file_path_glob))
    if len(tiff_paths) == 0:
        raise Exception("No files found at {}".format(file_path_glob))

    z_dim = len(tiff_paths)
    y_dim, x_dim = tifffile.imread(tiff_paths[0]).shape[0:2]

    return x_dim, y_dim, z_dim, tiff_paths


def crop_tiff_stack(
        tiff_paths: List[str],
        x_dim: int,
        y_dim: int,
        z_dim: int,
        output_filepath: str,
        output_tiff_stack: bool = False,
        crop_x: str = None,
        crop_y: str = None,
        crop_z: str = None) -> bool:
    """Crops TIFF stacks without loading the entire stack into memory at once
    Note: assumes at least one TIFF in the stack fits in memory

    Parameters
    ----------
    tiff_paths : List[str]
        A list of file path to the tiffs in the stack, in order of ascending z-coordinate
    x_dim : int
        The number of pixels or voxels in the x-dimension
    y_dim : int
        The number of pixels or voxels in the y-dimension
    z_dim : int
        The number of pixels or voxels in the z-dimension
    output_filepath : str
        A file path or directory to which the resulting TIFF file will be saved. If the
        file or directory does not exist, it will be created.
    output_tiff_stack : bool
        True if the cropped tiff should be saved as a collection of 2D tiffs,
        False if it should be saved in a multipage tiff
    crop_x : str
        A string in the format "[x_min],[x_max]" such that the resulting TIFF
        file is cropped to include all x-dimension pixels in the range
        (x_min, x_max), inclusive, of the original TIFF file
    crop_y : str
        A string in the format "[y_min],[y_max]" such that the resulting TIFF
        file is cropped to include all y-dimension pixels in the range
        (y_min, y_max), inclusive, of the original TIFF file
    crop_z : str
        A string in the format "[z_min],[z_max]" such that the resulting TIFF
        file is cropped to include all z-dimension pixels in the range
        (z_min, z_max), inclusive, of the original TIFF file

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if crop_x is None:
        x_min = 1
        x_max = x_dim
    else:
        x_min, x_max = [int(s) for s in crop_x.split(",")]
        if x_min < 0:
            raise ValueError("Minimum x-cropping value must be positive")
        if x_min > x_max:
            raise ValueError("Minimum x-cropping value must be smaller"
                             "than maximum x-cropping value")
        if x_max > x_dim:
            raise ValueError("Maximum x-cropping value must be smaller than"
                             "the dimensions of the image")

    if crop_y is None:
        y_min = 1
        y_max = y_dim
    else:
        y_min, y_max = [int(s) for s in crop_y.split(",")]
        if y_min < 0:
            raise ValueError("Minimum y-cropping value must be positive")
        if y_min > y_max:
            raise ValueError("Minimum y-cropping value must be smaller"
                             "than maxmimum y-cropping value")
        if y_max > y_dim:
            raise ValueError("Maximum y-cropping value must be smaller than"
                             "the dimensions of the image")

    if crop_z is None:
        z_min = 1
        z_max = z_dim
    else:
        z_min, z_max = [int(s) for s in crop_z.split(",")]
        if z_min < 0:
            raise ValueError("Minimum z-cropping value must be positive")
        if z_min > z_max:
            raise ValueError("Minimum z-cropping value must be smaller"
                             "than maxmimum z-cropping value")
        if z_max > z_dim:
            raise ValueError("Maximum z-cropping value must be smaller than"
                             "the dimensions of the image")

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

    return True


if __name__ == "__main__":
    main()
