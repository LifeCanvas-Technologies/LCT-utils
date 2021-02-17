"""
A module to enable command-line cropping of coordinates corresponding to tiffs, stored in json format
"""
# TODO: unit tests
import argparse
import sys
import json
from typing import List, Tuple, Union

try:
    from .utils import get_tiff_dims, parse_crop_dims_str, pretty_print_json_coords
except ModuleNotFoundError:
    from LCT_utils.utils import get_tiff_dims, parse_crop_dims_str, pretty_print_json_coords, validate_tiff_crop_dims


def parse_args(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Crop json coordinates")
    parser.add_argument("--input",
                        required=True,
                        help="A json file to crop")
    parser.add_argument("--output",
                        help="The path to the output json file (if not specified"
                             "will return number of coordinate tuples in input file)")
    parser.add_argument("--crop-x",
                        help="Specify as \"--crop-x [x-min],[x-max]\".")
    parser.add_argument("--crop-y",
                        help="Specify as \"--crop-y [y-min],[y-max]\".")
    parser.add_argument("--crop-z",
                        help="Specify as \"--crop-z [z-min],[z-max]\".")
    parser.add_argument("--flip-x",
                        help="Reverse the x-direction after scaling and cropping",
                        action="store_true")
    parser.add_argument("--flip-y",
                        help="Reverse the y-direction after scaling and cropping",
                        action="store_true")
    parser.add_argument("--flip-z",
                        help="Reverse the z-direction after scaling and cropping",
                        action="store_true")
    group_x_scale = parser.add_mutually_exclusive_group()
    group_x_scale.add_argument("--scale-x",
                               help="Scale the x-dimension by a given factor",
                               type=int,
                               default=1)
    group_x_scale.add_argument("--x-pixels",
                               help="The number of x-dimension pixels in the transformed tiff",
                               type=int)

    group_y_scale = parser.add_mutually_exclusive_group()
    group_y_scale.add_argument("--scale-y",
                               help="Scale the y-dimension by a given factor",
                               type=int,
                               default=1)
    group_y_scale.add_argument("--y-pixels",
                               help="The number of y-dimension pixels in the transformed tiff",
                               type=int)

    group_z_scale = parser.add_mutually_exclusive_group()
    group_z_scale.add_argument("--scale-z",
                               help="Scale the z-dimension by a given factor",
                               type=int,
                               default=1)
    group_z_scale.add_argument("--z-pixels",
                               help="The number of z-dimension pixels in the transformed tiff",
                               type=int)

    parser.add_argument("--original-tiff",
                        help="A UNIX-style glob to the original tiff or tiffs"
                             "required if x, y, or z flipping OR pixel count-based scaling is specified"
                             "AND the cropping dimensions are not specified.")

    parsed_args = parser.parse_args(args)

    if not parsed_args.original_tiff:
        if (parsed_args.flip_x and not parsed_args.crop_x) \
                or (parsed_args.flip_y and not parsed_args.crop_y)\
                or (parsed_args.flip_z and not parsed_args.crop_z):
            parser.error("--original-tiff is required when flipping is set on an axis for which cropping is not set")
        if (parsed_args.x_pixels and not parsed_args.crop_x) \
                or (parsed_args.y_pixels and not parsed_args.crop_y) \
                or (parsed_args.z_pixels and not parsed_args.crop_z):
            parser.error("--original-tiff is required when scaling based on pixel counting is set"
                         "for an axis for which cropping is not set")

    return parsed_args


def main(args=sys.argv[1:]):
    args = parse_args(args)
    if args.output is None:
        num_tuples = get_num_coords(input_path=args.input)
        print("JSON file at {} has {} coordinate tuples".format(args.input, num_tuples))
    else:
        if args.x_pixels or args.y_pixels or args.z_pixels:
            x_dim, y_dim, z_dim = get_tiff_dims(args.original_tiff, print_out=False)
            crop_xyz = parse_crop_dims_str(args.crop_x, args.crop_y, args.crop_z)
            crop_x, crop_y, crop_z = validate_tiff_crop_dims(x_dim, y_dim, z_dim, *crop_xyz)

            if args.x_pixels:
                args.scale_x = float(args.x_pixels) / (crop_x[1] - crop_x[0] + 1)
            if args.y_pixels:
                args.scale_y = float(args.y_pixels) / (crop_y[1] - crop_y[0] + 1)
            if args.z_pixels:
                args.scale_z = float(args.z_pixels) / (crop_z[1] - crop_z[0] + 1)

        transform_json(input_path=args.input,
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
                       original_tiff=args.original_tiff)


def transform_json(input_path: str,
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
                   original_tiff: str = None) -> None:
    """Writes json file of cropped coordinates given an input json and cropping specifications
    Assumes json field with coordinates is only content in input json file, maps old coordinates to
    new coordinates assuming that the lower bounds on each of crop_x, crop_y, and crop_z will be (0, 0, 0)

    Also assumes coordinates are 3D

    Order of operations: cropping > scaling > flipping

    Parameters
    ----------
    input_path : str
        Path to input json file with uncropped coordinates
    output_path : str
        Path to output the resulting json file
    crop_x : Union[str, List[int], Tuple[int, int]]
        String in format "[x_min,x_max]" or Sequence of ints in format [x_min, x_max] specifying the minimum
        and maximum x-values (inclusive) to include in the resulting coordinates
    crop_y : Union[str, List[int], Tuple[int, int]]
        String in format "[y_min,y_max]" or Sequence of ints in format [y_min, y_max] specifying the minimum
        and maximum y-values (inclusive) to include in the resulting coordinates
    crop_z : Union[str, List[int], Tuple[int, int]]
        String in format "[z_min,z_max]" or Sequence of ints in format [z_min, z_max] specifying the minimum
        and maximum x-values (inclusive) to include in the resulting coordinates
    flip_x : bool
        Specify whether or not the output coordinates should be reversed along the x-axis after cropping
    flip_y : bool
        Specify whether or not the output coordinates should be reversed along the y-axis after cropping
    flip_z : bool
        Specify whether or not the output coordinates should be reversed along the z-axis after cropping
    scale_x : float
        Specify the scaling factor for x-coordinates after cropping
    scale_y : float
        Specify the scaling factor for y-coordinates after cropping
    scale_z : float
        Specify the scaling factor for z-coordinates after cropping
    original_tiff : str
        UNIX-style glob expression to original TIFF stack
        Required for flipping

    Returns
    -------
    None
    """
    crop_x, crop_y, crop_z = parse_crop_dims_str(crop_x=crop_x, crop_y=crop_y, crop_z=crop_z)
    if (crop_x[1] == float("Inf") and flip_x) or \
            (crop_y[1] == float("Inf") and flip_y) or \
            (crop_z[1] == float("Inf") and flip_z):
        x_dim, y_dim, z_dim = get_tiff_dims(original_tiff)
        if flip_x:
            crop_x = 1, x_dim
        if flip_y:
            crop_y = 1, y_dim
        if flip_z:
            crop_z = 1, z_dim

    with open(input_path) as fd:
        full_coords = json.load(fd)
    print("Reading in {} coordinate tuples from {}".format(len(full_coords), input_path))

    cropped_coords = []
    for x, y, z in full_coords:
        if crop_x[0] <= x <= crop_x[1] and crop_y[0] <= y <= crop_y[1] and crop_z[0] <= z <= crop_z[1]:

            # Cropping
            if crop_x[0] != -float("Inf"):
                cropped_x = x - crop_x[0]
            else:
                cropped_x = x
            if crop_y[0] != -float("Inf"):
                cropped_y = y - crop_y[0]
            else:
                cropped_y = y
            if crop_z[0] != -float("Inf"):
                cropped_z = z - crop_z[0]
            else:
                cropped_z = z

            # Scaling
            cropped_x = scale_x * cropped_x
            cropped_y = scale_y * cropped_y
            cropped_z = scale_z * cropped_z

            # Flipping
            if flip_x or flip_y or flip_z:
                if original_tiff is None:
                    raise ValueError
                else:
                    x_dim, y_dim, z_dim = get_tiff_dims(original_tiff)
                    if flip_x:
                        cropped_x = x_dim - cropped_x + 1
                    if flip_y:
                        cropped_y = y_dim - cropped_y + 1
                    if flip_z:
                        cropped_z = z_dim - cropped_z + 1

            cropped_coord = [cropped_x, cropped_y, cropped_z]
            cropped_coords.append(cropped_coord)

    print("Writing {} coordinate tuples to {}".format(len(cropped_coords), output_path))
    pretty_print_json_coords(cropped_coords, output_path)


def get_num_coords(input_path: str) -> int:
    """Simple function to return the number of coordinate tuples from a given json file

    Parameters
    ----------
    input_path : str
        Path to the json file

    Returns
    -------
    int
    """
    with open(input_path) as fd:
        return len(json.load(fd))


if __name__ == "__main__":
    main()
