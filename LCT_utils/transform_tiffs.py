"""
A module to enable command-line cropping of TIFF files and stacks
"""
# TODO: rotate image 90 degrees
# ^ order with cropping/flipping/scaling?
# Trivial to do with single tiff / roll of TIFF stack--pitch/yaw seems incredibly difficult with TIFF stack
# TODO: unit tests
# TODO: specify n-cores to use

import argparse
import sys
import os
import platform
from typing import List, Tuple, Sequence, Union
import multiprocessing as mp
import warnings

import numpy as np
from scipy.ndimage import zoom
from skimage import img_as_uint
import tifffile
import tqdm.auto as tqdm


try:
    from .utils import \
        glob_to_list, get_tiff_dims, parse_crop_dims_str, validate_tiff_crop_dims, reshape_single_tiff
except ModuleNotFoundError:
    from LCT_utils.utils import \
        glob_to_list, get_tiff_dims, parse_crop_dims_str, validate_tiff_crop_dims, reshape_single_tiff


def parse_args(args=sys.argv[1:]):
    """Reads in command line arguments to script"""
    parser = argparse.ArgumentParser(description="Transform tiff files.")
    parser.add_argument("--input",
                        required=True,
                        help="A tiff file to crop.")
    parser.add_argument("--output",
                        help="The path to the output tiff file (if not "
                             "specified, will return dimensions of input tiff).")
    parser.add_argument("--crop-x",
                        help="The minimum and maximum x-values in the original tiff to include in the transformed tiff."
                             "Specify as \"--crop-x [x-min],[x-max]\".")
    parser.add_argument("--crop-y",
                        help="The minimum and maximum y-values in the original tiff to include in the transformed tiff."
                             "Specify as \"--crop-y [y-min],[y-max]\".")
    parser.add_argument("--crop-z",
                        help="The minimum and maximum z-values in the original tiff to include in the transformed tiff."
                             "Specify as \"--crop-z [z-min],[z-max]\".")
    parser.add_argument("--flip-x",
                        help="Reverse the x-direction of the image after scaling and cropping",
                        action="store_true")
    parser.add_argument("--flip-y",
                        help="Reverse the y-direction of the image after scaling and cropping",
                        action="store_true")
    parser.add_argument("--flip-z",
                        help="Reverse the z-direction of the image after scaling and cropping",
                        action="store_true")
    group_x_scale = parser.add_mutually_exclusive_group()
    group_x_scale.add_argument("--scale-x",
                               help="Scale the x-dimension by a given factor",
                               type=float,
                               default=1)
    group_x_scale.add_argument("--x-pixels",
                               help="The number of x-dimension pixels in the transformed tiff",
                               type=int)

    group_y_scale = parser.add_mutually_exclusive_group()
    group_y_scale.add_argument("--scale-y",
                               help="Scale the y-dimension by a given factor",
                               type=float,
                               default=1)
    group_y_scale.add_argument("--y-pixels",
                               help="The number of y-dimension pixels in the transformed tiff",
                               type=int)

    group_z_scale = parser.add_mutually_exclusive_group()
    group_z_scale.add_argument("--scale-z",
                               help="Scale the z-dimension by a given factor",
                               type=float,
                               default=1)
    group_z_scale.add_argument("--z-pixels",
                               help="The number of z-dimension pixels in the transformed tiff",
                               type=int)
    parser.add_argument("--split-output",
                        help="Specify whether or not the output should be saved as a collection of tiffs"
                             "Note: the OUTPUT should be a path to a folder, not a file.",
                        action="store_true")
    parser.add_argument("--start-num",
                        help="Specify the number at which to start split output TIFFs on.",
                        type=int,
                        default=0)
    parser.add_argument("--num-digits",
                        help="Specify the number of digits in the output files",
                        type=int)

    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=sys.argv[1:]):
    args = parse_args(args)
    if args.output is None:
        get_tiff_dims(input_glob_or_list=args.input,
                      print_out=True)
    else:
        if args.x_pixels or args.y_pixels or args.z_pixels:
            x_dim, y_dim, z_dim = get_tiff_dims(args.input, print_out=False)
            crop_xyz = parse_crop_dims_str(args.crop_x, args.crop_y, args.crop_z)
            crop_x, crop_y, crop_z = validate_tiff_crop_dims(x_dim, y_dim, z_dim, *crop_xyz)

            if args.x_pixels:
                args.scale_x = float(args.x_pixels) / (crop_x[1] - crop_x[0] + 1)
            if args.y_pixels:
                args.scale_y = float(args.y_pixels) / (crop_y[1] - crop_y[0] + 1)
            if args.z_pixels:
                args.scale_z = float(args.z_pixels) / (crop_z[1] - crop_z[0] + 1)

        transform_tiff(input_glob_or_list=args.input,
                       output_path=args.output,
                       crop_x=args.crop_x,
                       crop_y=args.crop_y,
                       crop_z=args.crop_z,
                       flip_x=args.flip_x,
                       flip_y=args.flip_y,
                       flip_z=args.flip_z,
                       scale_x=args.scale_x,
                       scale_y=args.scale_y,
                       scale_z=args.scale_z,
                       split_output=args.split_output,
                       start_num=args.start_num,
                       num_digits=args.num_digits)


def transform_tiff(input_glob_or_list: Union[str, List[str]],
                   output_path: str,
                   crop_x: Union[str, List[int], Tuple[int, int]] = None,
                   crop_y: Union[str, List[int], Tuple[int, int]] = None,
                   crop_z: Union[str, List[int], Tuple[int, int]] = None,
                   flip_x: bool = False,
                   flip_y: bool = False,
                   flip_z: bool = False,
                   scale_x: float = 1,
                   scale_y: float = 1,
                   scale_z: float = 1,
                   split_output: bool = False,
                   start_num: int = 0,
                   num_digits: int = None) -> None:
    """Main function for cropping tiffs

    Note: progress bars not particularly accurate for inputs with very large tiff files due to
    IO-limited multiprocessing and how the jobs are queued

    Parameters
    ----------
    input_glob_or_list : Union[str, List[str]]
        String containing input glob or list of file paths to individual tiffs
    output_path : str
        Output file path, should either be a directory or tiff file, depending on the value of split_output
    crop_x : Union[str, List[int], Tuple[int, int]]
        String in format "[x_min,x_max]" or Sequence of ints in format [x_min, x_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff
    crop_y : Union[str, List[int], Tuple[int, int]]
        String in format "[y_min,y_max]" or Sequence of ints in format [y_min, y_max] specifying the minimum
        and maximum y-values (inclusive) to include in the outputted cropped tiff
    crop_z : Union[str, List[int], Tuple[int, int]]
        String in format "[z_min,z_max]" or Sequence of ints in format [z_min, z_max] specifying the minimum
        and maximum x-values (inclusive) to include in the outputted cropped tiff
    flip_x : bool
        Specify whether or not the output image should have its x-axis reversed after cropping
    flip_y : bool
        Specify whether or not the output image should have its y-axis reversed after cropping
    flip_z : bool
        Specify whether or not the output image should have its z-axis reversed after cropping
    scale_x : float
        Specify the scaling factor along the x-axis after cropping
    scale_y : float
        Specify the scaling factor along the y-axis after cropping
    scale_z : float
        Specify the scaling factor along the z-axis after cropping
    split_output : bool
        True if output should be a directory of 2D tiffs, False if output should be single multi-page tiff file
    start_num : int
        Number from which to start outputting in output TIFF stack
    num_digits : int
        Number of digits in output TIFF stack file names

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

    x_dim, y_dim, z_dim = get_tiff_dims(input_glob_or_list)
    crop_x, crop_y, crop_z = parse_crop_dims_str(crop_x, crop_y, crop_z)
    crop_x, crop_y, crop_z = validate_tiff_crop_dims(x_dim, y_dim, z_dim, crop_x, crop_y, crop_z)

    if num_digits is None:
        num_digits = int(np.floor(np.log10(crop_z[1] - crop_z[0] + 1)) + 1)

    if len(tiff_paths) == 1:
        tiff_ndarray_reshape = reshape_single_tiff(tiff_paths[0])
        transform_single_tiff(tiff_ndarray_reshape,
                              output_filepath=output_path,
                              output_tiff_stack=split_output,
                              crop_x=crop_x,
                              crop_y=crop_y,
                              crop_z=crop_z,
                              flip_x=flip_x,
                              flip_y=flip_y,
                              flip_z=flip_z,
                              scale_x=scale_x,
                              scale_y=scale_y,
                              scale_z=scale_z,
                              start_num=start_num,
                              num_digits=num_digits)
    else:
        transform_multi_tiff(tiff_paths,
                             output_filepath=output_path,
                             output_tiff_stack=split_output,
                             crop_x=crop_x,
                             crop_y=crop_y,
                             crop_z=crop_z,
                             flip_x=flip_x,
                             flip_y=flip_y,
                             flip_z=flip_z,
                             scale_x=scale_x,
                             scale_y=scale_y,
                             scale_z=scale_z,
                             start_num=start_num,
                             num_digits=num_digits)


def transform_single_tiff(
        tiff_ndarray: np.ndarray,
        output_filepath: str,
        output_tiff_stack: bool,
        crop_x: Tuple[int, int],
        crop_y: Tuple[int, int],
        crop_z: Tuple[int, int],
        flip_x: bool,
        flip_y: bool,
        flip_z: bool,
        scale_x: float,
        scale_y: float,
        scale_z: float,
        num_digits: int,
        start_num: int = 0) -> None:
    """Crops a TIFF ndarray, saving as a TIFF file

    Order of operations: cropping > scaling > flipping

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
        False if it should be saved in a multi-page tiff
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
    flip_x : bool
        Specify whether or not the output image should have its x-axis reversed after cropping
    flip_y : bool
        Specify whether or not the output image should have its y-axis reversed after cropping
    flip_z : bool
        Specify whether or not the output image should have its z-axis reversed after cropping
    scale_x : float
        Specify the scaling factor along the x-axis after cropping
    scale_y : float
        Specify the scaling factor along the y-axis after cropping
    scale_z : float
        Specify the scaling factor along the z-axis after cropping
    start_num : int
        Number from which to start outputting in output TIFF stack
    num_digits : int
        Number of digits in output TIFF stack file names

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
    tiff_ndarray_scaled = zoom(tiff_ndarray_cropped, (scale_z, scale_y, scale_x))
    if flip_z:
        tiff_ndarray_scaled = np.flip(tiff_ndarray_scaled, 0)
    if flip_y:
        tiff_ndarray_scaled = np.flip(tiff_ndarray_scaled, 1)
    if flip_x:
        tiff_ndarray_scaled = np.flip(tiff_ndarray_scaled, 2)

    if output_tiff_stack is False:
        with tifffile.TiffWriter(output_filepath, bigtiff=bigtiff) as file:
            file.save(tiff_ndarray_scaled)
    else:
        if not os.path.exists(output_filepath):
            print("Creating directory at {}.".format(output_filepath))
            os.mkdir(output_filepath)
        for tiff_num, z in enumerate(tqdm.tqdm(range(tiff_ndarray_scaled.shape[0]))):
            filename = output_filepath + \
                       "img_" + \
                       str(tiff_num + start_num).zfill(num_digits) \
                       + ".tiff"
            with tifffile.TiffWriter(filename, bigtiff=bigtiff) as file:
                file.save(tiff_ndarray_cropped[z])


def transform_multi_tiff(
        tiff_paths: List[str],
        output_filepath: str,
        output_tiff_stack: bool,
        crop_x: Sequence[int],
        crop_y: Sequence[int],
        crop_z: Sequence[int],
        flip_x: bool,
        flip_y: bool,
        flip_z: bool,
        scale_x: float,
        scale_y: float,
        scale_z: float,
        num_digits: int,
        start_num: int = 0) -> None:
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
    flip_x : bool
        Specify whether or not the output image should have its x-axis reversed after cropping
    flip_y : bool
        Specify whether or not the output image should have its y-axis reversed after cropping
    flip_z : bool
        Specify whether or not the output image should have its z-axis reversed after cropping
    scale_x : float
        Specify the scaling factor along the x-axis after cropping
    scale_y : float
        Specify the scaling factor along the y-axis after cropping
    scale_z : float
        Specify the scaling factor along the z-axis after cropping
    start_num : int
        Number from which to start outputting in output TIFF stack
    num_digits : int
        Number of digits in output TIFF stack file names

    Returns
    -------
    None
    """
    x_min, x_max = crop_x
    y_min, y_max = crop_y
    z_min, z_max = crop_z

    z_pixels = scale_z * (z_max - z_min + 1)
    z_range = np.linspace(z_min - 1, z_max - 1, np.round(z_pixels))
    if flip_z:
        z_range = np.flip(z_range)

    if not output_tiff_stack:
        sys.exit("DEV: not supported currently")

    if not os.path.exists(output_filepath):
        print("Creating directory at {}.".format(output_filepath))
        os.mkdir(output_filepath)

    with mp.Pool() as pool:
        results = []
        for tiff_num, z in enumerate(z_range):
            filename = output_filepath + "img_" + \
                       str(tiff_num + start_num).zfill(num_digits) + ".tiff"
            if z.is_integer():
                result = pool.apply_async(write_one_plane, ([tiff_paths[int(z)]],
                                                            filename,
                                                            x_min,
                                                            x_max,
                                                            y_min,
                                                            y_max,
                                                            scale_x,
                                                            scale_y,
                                                            flip_x,
                                                            flip_y))
            else:
                lower_bound, upper_bound = int(np.floor(z)), int(np.ceil(z))
                result = pool.apply_async(write_two_planes, (tiff_paths[lower_bound:upper_bound + 1],
                                                             filename,
                                                             z - lower_bound,
                                                             x_min,
                                                             x_max,
                                                             y_min,
                                                             y_max,
                                                             scale_x,
                                                             scale_y,
                                                             flip_x,
                                                             flip_y))
            results.append(result)

        [result.get() for result in tqdm.tqdm(results)]


def write_one_plane(filepath: str,
                    outpath: str,
                    x_min: int,
                    x_max: int,
                    y_min: int,
                    y_max: int,
                    scale_x: float,
                    scale_y: float,
                    flip_x: bool,
                    flip_y: bool):
    """Writes out a cropped, scaled, and flipped single TIFF file from single TIFF
    Order of operations: cropping, scaling, flipping

    Parameters
    ----------
    filepath : str
        File path of original TIFF
    outpath : str
        File path to output TIFF
    x_min : int
        Minimum x-coordinate to keep in the output TIFF file
    x_max : int
        Maximum x-coordinate to keep in the output TIFF file
    y_min : int
        Minimum y-coordinate to keep in the output TIFF file
    y_max : int
        Maximum y-coordinate to keep in the output TIFF file
    scale_x : float
        Scaling factor in the x-dimension
    scale_y : float
        Scaling factor in the y-dimension
    flip_x : bool
        True if the output image should be reflected about the y axis (reversing the image along the x-dimension)
    flip_y : bool
        True if the output image should be reflected about the x axis (reversing the image along the y-dimension)

    Returns
    -------
    None
    """
    plane_ndarray = img_as_uint(tifffile.imread(filepath)[y_min - 1:y_max, x_min - 1:x_max])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plane = img_as_uint(zoom(plane_ndarray, (scale_y, scale_x)))

    if flip_x:
        plane = np.flip(plane, 1)
    if flip_y:
        plane = np.flip(plane, 0)

    with tifffile.TiffWriter(outpath) as file:
        file.save(plane)


def write_two_planes(filepaths: List[str],
                     outpath: str,
                     z_level: float,
                     x_min: int,
                     x_max: int,
                     y_min: int,
                     y_max: int,
                     scale_x: float,
                     scale_y: float,
                     flip_x: bool,
                     flip_y: bool):
    """Writes out a cropped, scaled, and flipped TIFF file interpolated between two TIFFs
    Order of operations: cropping, scaling, flipping

    Parameters
    ----------
    filepaths : List[str]
        List of file paths of original two TIFFs
    outpath : str
        File path to output TIFF
    z_level : float
        Float between 0 and 1 describing interpolation position between two planes, 0 being close to bottom
    x_min : int
        Minimum x-coordinate to keep in the output TIFF file
    x_max : int
        Maximum x-coordinate to keep in the output TIFF file
    y_min : int
        Minimum y-coordinate to keep in the output TIFF file
    y_max : int
        Maximum y-coordinate to keep in the output TIFF file
    scale_x : float
        Scaling factor in the x-dimension
    scale_y : float
        Scaling factor in the y-dimension
    flip_x : bool
        True if the output image should be reflected about the y axis (reversing the image along the x-dimension)
    flip_y : bool
        True if the output image should be reflected about the x axis (reversing the image along the y-dimension)

    Returns
    -------
    None
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p1 = img_as_uint(zoom(tifffile.imread(filepaths[0])[y_min - 1:y_max, x_min - 1:x_max], (scale_y, scale_x)))
        p2 = img_as_uint(zoom(tifffile.imread(filepaths[1])[y_min - 1:y_max, x_min - 1:x_max], (scale_y, scale_x)))

    # Uses linear spline to interpolate between two planes
    p_comb = ((1 - z_level) * p1) + (z_level * p2)
    p_comb = img_as_uint(np.around(p_comb).astype(np.uint16))

    if flip_x:
        p_comb = np.flip(p_comb, 1)
    if flip_y:
        p_comb = np.flip(p_comb, 0)

    with tifffile.TiffWriter(outpath) as file:
        file.save(p_comb)


if __name__ == "__main__":
    main()
