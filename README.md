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
### The transform-tiff command
**transform-tiff** is the main function of LCT_utils. It reads in a specified TIFF file (2D or 3D) or 
UNIX-style glob expression pointing to a collection of 2D TIFF files and outputs a TIFF file 
or collection of TIFF files cropped, reflected, and scaled to specification. 
    
**transform-tiff** is run as follows:
    
    transform-tiff 
        --input <input TIFF path or glob> \
        [--output <output TIFF path or directory>] \
        [--crop-x <x-min>, <x-max>] \
        [--crop-y <y-min>, <y-max>] \
        [--crop-z <z-min>, <z-max>] \
        [--flip-x] \
        [--flip-y] \
        [--flip-z] \
        [--scale-x <scaling factor> | --x-pixels <num pixels>] \
        [--scale-y <scaling factor> | --y-pixels <num pixels>] \
        [--scale-z <scaling factor> | --z-pixels <num pixels>] \        
        [--split-output]

where:
* *--input* is a path to a TIFF file or a UNIX-style glob expression pointing to a collection of TIFFs 
* *--output* is the path to which the cropped TIFF file (or a directory if *--split-output* is included)
 should be written; if *--output* is not included, **crop-tiff** will instead print out the dimensions
 of the input TIFF file(s)
* *--crop-x*, *--crop-y*, and *--crop-z* specify the minimum and maximum values in the 
x-, y-, and z-dimensions, respectively, that should be included in the crop
* *--flip-x*, *--flip-y*, and *--flip-z* can be included to reverse the image in the 
x-, y-, and z-directions, respectively
* *--scale-x*, *--scale-y*, and *--scale-z* specify a scaling factor for the 
output tiff file in the x-, y-, and z-direction, respectively
* *--x-pixels*, *--y-pixels*, *--z-pixels* can be used in place of the scaling factors to 
specify the number of pixels in the x-, y-, and z-dimensions, respectively, in the output tiff
* *--split-output* specifies that the resulting cropped TIFF file should be written out as stack of
2D TIFFs at a directory specified by *--output*. These TIFF files will automatically be numbered
in the format "img_[####].tiff", where [####] specifies an automatically generated sequential number,
specifying that TIFF's relative z-position in the stack. 


### The transform-json command
**transform-json** reads in a specified JSON file of 3D coordinates and trims them such that they 
remain accurate for a TIFF file cropped to the same specifications

**transform-json** is used as follows:

    transform-json
        --input <input JSON path> \
        [--output <output JSON path>] \
        [--crop-x <x-min>, <x-max>] \
        [--crop-y <y-min>, <y-max>] \
        [--crop-z <z-min>, <z-max>] \
        [--flip-x] \
        [--flip-y] \
        [--flip-z] \
        [--scale-x <scaling factor> | --x-pixels <num pixels>] \
        [--scale-y <scaling factor> | --y-pixels <num pixels>] \
        [--scale-z <scaling factor> | --z-pixels <num pixels>]

where:
* *--input* is a path to a JSON file of 3D coordinates **only** 
* *--output* is the path to which the cropped JSON file should be written; 
if *--output* is not included, **crop-json** will instead print out the number of coordinate triples
are in the input JSON file
* *--crop-x*, *--crop-y*, and *--crop-z* specify the minimum and maximum values in the 
x-, y-, and z-dimensions, respectively, that should be included in the crop
* *--flip-x*, *--flip-y*, and *--flip-z* can be included to reverse the coordinates along the 
x-, y-, and z-directions, respectively
* *--scale-x*, *--scale-y*, and *--scale-z* specify a scaling factor for the 
output coordinates in the x-, y-, and z-direction, respectively
* *--x-pixels*, *--y-pixels*, *--z-pixels* can be used in place of the scaling factors to 
specify the number of pixels in the x-, y-, and z-dimensions, respectively, in the corresponding 
output tiff


### The json-to-tiff command
**json-to-tiff** reads in a specified JSON file of 3D coordinates and turns them into an Imaris 
compatible TIFF stack for visualization. 

**json-to-tiff** is used as follows:

    json-to-tiff
        --input <input JSON path> \
        --output <output TIFF directory> \
        --reference-tiff <original TIFF file or stack>

where:
* *--input* is a path to a JSON file of 3D coordinates
* *--output* is a directory to where the output TIFF stack should be written
* *--reference-tiff* is a UNIX-style glob expression to the original TIFF stack from which the
input JSON coordinates are derived.  


## Import Usage
**LCT-utils** can also be used as a Python package through import statements. To do so, we recommend:
   
    import LCT_utils as lct
 
Then, `lct.transform_tiff()`, `lct.get_tiff_dims()`, 
`lct.transform_json()`, `lct.get_num_coords()` are imported as functions with similar usages to above. 
See the source code for detailed documentation. 
