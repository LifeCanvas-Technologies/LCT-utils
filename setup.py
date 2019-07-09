from setuptools import setup

version = "0.1.0"
description = "A Python package that enables command-line cropping of TIFF files and stacks"

setup(
    name="crop_tiffs",
    version=version,
    description=description,
    install_requires=[
        "numpy",
        "tifffile",
        'tqdm'
    ],
    author="Richard Qiu",
    packages=["crop_tiffs"],
    entry_points={
        "console_scripts": ["crop-tiff=crop_tiffs.crop_tiffs:main"]
    }
)
