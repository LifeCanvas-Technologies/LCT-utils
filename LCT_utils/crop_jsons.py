"""
A module to enable command-line cropping of coordinates corresponding to tiffs, stored in json format
"""

import argparse
import sys
import json
from typing import List, Tuple, Union
from .utils import parse_crop_dims_str, pretty_print_json_coords


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

    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=sys.argv[1:]):
    args = parse_args(args)
    if args.output is None:
        num_tuples = get_num_tuples(input_path=args.input)
        print("JSON file at {} has {} coordinate tuples".format(args.input, num_tuples))
    else:
        crop_json(input_path=args.input,
                  output_path=args.output,
                  crop_x=args.crop_x,
                  crop_y=args.crop_y,
                  crop_z=args.crop_z)


def crop_json(input_path: str,
              output_path: str,
              crop_x: Union[str, List[int], Tuple[int, int]] = None,
              crop_y: Union[str, List[int], Tuple[int, int]] = None,
              crop_z: Union[str, List[int], Tuple[int, int]] = None) -> None:
    """Writes json file of cropped coordinates given an input json and cropping specifications
    Assumes json field with coordinates is only content in input json file, maps old coordinates to
    new coordinates assuming that the lower bounds on each of crop_x, crop_y, and crop_z will be (0, 0, 0)

    Also assumes coordinates are 3D

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

    Returns
    -------
    None
    """
    crop_x, crop_y, crop_z = parse_crop_dims_str(crop_x=crop_x, crop_y=crop_y, crop_z=crop_z)

    with open(input_path) as fd:
        full_coords = json.load(fd)
    print("Reading in {} coordinate tuples from {}".format(len(full_coords), input_path))

    cropped_coords = []
    for x, y, z in full_coords:
        if crop_x[0] <= x <= crop_x[1] and crop_y[0] <= y <= crop_y[1] and crop_z[0] <= z <= crop_z[1]:
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
            cropped_coord = [cropped_x, cropped_y, cropped_z]
            cropped_coords.append(cropped_coord)

    print("Writing {} coordinate tuples to {}".format(len(cropped_coords), output_path))
    pretty_print_json_coords(cropped_coords, output_path)


def get_num_tuples(input_path: str) -> int:
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
