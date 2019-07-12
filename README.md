# LCT_utils
LCT_utils is a small utility package that enables the command-line cropping of very large TIFF files
and particularly formatted JSON files with coordinates corresponding to TIFF files. 

In particular, this package was predominantly designed to ease the handling of large quantities of
volumetric data produced in the process of high-throughput microscopy, but the modules making up LCT_utils
can easily be adapted to other use cases. 

## Installation
LCT_utils can be install using pip:

`$ pip install LCT_utils`

Alternatively, a dev version of LCT_utils can be installed by running:

`$ pip install git+https://github.com/lifecanvastech/LCT-utils.git`

## Command Line Usage
### The crop-tiff command
**crop-tiff** is the main function of LCT_utils. It reads in a specified TIFF file (2D or 3D) or UNIX-style
glob expression pointing to a collection of 2D TIFF files. 
    
**crop-tiff** is run as follows:
    
    crop-tiff 
        --input <input TIFF path or glob> \
        [--output <output TIFF path or directory>] \
        [--crop-x <x-min>, <x-max>] \
        [--crop-y <y-min>, <y-max>] \
        [--crop-z <z-min>, <z-max>] \
        [--split-output]

where:
* *--input* is a path to a TIFF file or a UNIX-style glob expression pointing to a collection of TIFFs 
* *--output* is the path to which the cropped TIFF file (or a directory if *--split-output* is included)
 should be written; if *--output* is not included, **crop-tiff** will instead print out the dimensions
 of the input TIFF file(s)
* *--crop-x* specifies the minimum and maximum values in the x-dimension that should be included
in the crop
* *--crop-y* specifies the minimum and maximum values in the y-dimension that should be included
in the crop
* *--crop-z* specifies the minimum and maximum values in the z-dimension that should be included
in the crop
* *--split-output* specifies that the resulting cropped TIFF file should be written out as stack of
2D TIFFs at a directory specified by *--output*. These TIFF files will automatically be numbered
in the format "img_[####].tiff", where [####] specifies an automatically generated sequential number,
specifying that TIFF's relative z-position in the stack. 


### The crop-json command
**crop-json** reads in a specified JSON file of 3D coordinates and trims them such that they 
remain accurate for a TIFF file cropped to the same specifications

**crop-json** is used as follows:

    crop-json
        --input <input JSON path> \
        [--output <output JSON path>] \
        [--crop-x <x-min>, <x-max>] \
        [--crop-y <y-min>, <y-max>] \
        [--crop-z <z-min>, <z-max>] \

where:
* *--input* is a path to a JSON file of 3D coordinates **only** 
* *--output* is the path to which the cropped JSON file should be written; 
if *--output* is not included, **crop-json** will instead print out the number of coordinate triples that
are in the input JSON file
* *--crop-x* specifies the minimum and maximum values in the x-dimension that should be included
in the crop
* *--crop-y* specifies the minimum and maximum values in the y-dimension that should be included
in the crop
* *--crop-z* specifies the minimum and maximum values in the z-dimension that should be included
in the crop


## Import Usage
**LCT-utils** can also be used as a Python package through import statements. To do so, we recommend:
   
    import LCT_utils as lct
 
 Then, `lct.crop_tiff()`, `lct.get_tiff_dims()`, `lct.crop_json()`, `lct.get_num_coords()` are imported
 as functions with similar usages to above. See the source code for detailed documentation. 
