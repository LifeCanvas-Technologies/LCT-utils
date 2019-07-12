# LCT_utils
LCT_utils is a small utility package that enables the command-line cropping of very large TIFF files
and particularly formatted JSON files with coordinates corresponding to TIFF files. 

In particular, this package was predominantly designed to ease the handling of large quantities of
volumetric data produced in the process of high-throughput microscopy, but the modules making up LCT_utils
can easily be adapted to other use cases. 

## Installation
LCT_utils can be install using pip:

`$ pip install LCT_utils`

Alternatively, LCT_utils can be installed by downloading this GitHub repo, `cd`'ing to the unzipped 
directory and running:

`$ pip install .`

## Command Line Usage
### The crop-tiff command
**crop-tiff** is the main function of LCT_utils. It reads in a specified TIFF file (2D or 3D) or UNIX-style
glob expression pointing to a collection of 2D TIFF files. 
    
**crop-tiff** is run as follows:
    
    crop-tiff 
        --input <input TIFF path or glob> \
        --input 

## Import Usage
