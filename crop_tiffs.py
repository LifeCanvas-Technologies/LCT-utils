"""
A module to enable command-line croppping of TIFF files and stacks
"""

import argparse
import glob
import sys
import numpy as n
import tifffile


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
	parser.add_argument("--large",
						help="Specify whether file to be cropped is too large to fit in memory"
							"Note: currently only works with TIFF stacks of 2D TIFFs.")

	parsed_args = parser.parse_args(args)
	return parsed_args


def main(args=sys.argv[1:]):
	args = _parse_args(args)

	try:
		tiff_ndarray = read_tiffs(args.input)
	except MemoryError as e:
		# TODO: write error to log file
		print("Error: not enough memory to load entire TIFF.")
		print("Use \"--large\" tag to attempt handling this TIFF.")
		exit(2)
	x_dim, y_dim, z_dim, tiff_ndarray_reshape = reshape_and_get_tiff_dims(
		tiff_ndarray)

	if args.output is None:
		if z_dim == 1:
			units_str = "[units: pixels]"
		else:
			units_str = "[units: voxels]"
		dims_str = str(x_dim) + ", " + str(y_dim) + ", " + str(z_dim)
		print("Dimensions (x, y, z) (width, height, depth)",
			units_str, " : ",  dims_str)
	else:
		crop_tiffs(tiff_ndarray_reshape,
				   args.output,
				   x_dim = x_dim,
				   y_dim = y_dim,
				   z_dim = z_dim,
				   crop_x = args.crop_x,
				   crop_y = args.crop_y,
				   crop_z = args.crop_z)


def read_tiffs(file_path_glob: str) -> numpy.ndarray:
	"""Reads in TIFF file from UNIX-style glob, returns numpy ndarray

	Parameters
	----------
	file_path_glob : str
		A glob expression for the stack of images to analyze,
		e.g. "/path-to-tiffs/*.tiff"

		If using a collection of 2D TIFFs, files must be in alphabetical
		order (e.g. by padding TIFF file names such as "image_00001.tiff")

	Returns
	-------
	numpy.ndarray
		Returns TIFF file(s) in corresponding numpy array format
		If multiple files or tiff stack, corresponding array dims are:
		<z, x, y>
		If extra channels included, will be returned as last array dim
	"""
	tiff_paths = sorted(glob.glob(file_path_glob))
	tiff_ndarray = tifffile.imread(tiff_paths)

	return tiff_ndarray


def reshape_and_get_tiff_dims(tiff_ndarray: numpy.ndarray) -> numpy.ndarray:
	"""Given a tiff ndarray, returns reformatted ndarray and its dimensions

	Note: assumes there are no more than 9 extra channels and that
	there are at least 10 pixels in the y dimension

	Parameters
	----------
	tiff_ndarray : numpy.ndarray
		A numpy ndarray representing a TIFF file or TIFF stack; extra
		channels such as alpha, RGB/CMYK, etc., are acceptable

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
		if tiff_ndarray.shape[-1] < 10:
			# Interpret as 3D TIFF with extra channels
			tiff_ndarray_reshape = tiff_ndarray.copy()
		else:
			sys.stderr.write("Not yet configured to handle this kind of TIFF.")
			exit(1)
	else:
		sys.stderr.write("Not yet configured to handle this kind of TIFF.")
		exit(1)

	z_dim, y_dim, x_dim = tiff_ndarray_reshape.shape[0:3]
	return x_dim, y_dim, z_dim, tiff_ndarray_reshape


def crop_tiffs(
		tiff_ndarray: numpy.ndarray,
		output_filepath: str,
		x_dim: int,
		y_dim: int,
		z_dim: int,
		crop_x=None: str,
		crop_y=None: str,
		crop_z=None: str) -> bool:
	"""Crops a TIFF ndarray, saving as a TIFF file
	Parameters
	----------
	tiff_ndarray : numpy.ndarray
		A numpy ndarray, already padded as needed to include 3 or 4
		dimensions, depending on whether or not there are extra channels,
		in the following order: [z, x, y, <extra channels>]
	output_filepath : str
		A file path to which the resulting TIFF file will be saved. If the
		file does not exist, it will be created.
	x_dim : int
		The number of pixels or voxels in the x-dimension
	y_dim : int
		The number of pixels or voxels in the y-dimension
	z_dim : int
		The number of pixels or voxels in the z-dimension
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
		x_min = 0
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
		y_min = 0
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
		z_min = 0
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

	tiff_ndarray_cropped = tiff_ndarray[z_min:z_max+1, \
										   y_min:y_max+1, \
										   x_min:x_max+1]

	with tifffile.TiffWriter(output_filepath) as file:
		file.save(tiff_ndarray_cropped)

	return True


if __name__ == "__main__":
	main()

