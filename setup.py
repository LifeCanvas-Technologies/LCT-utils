from setuptools import setup

version = "0.2.2"
description = "Command-line cropping of TIFF files and stacks as well as corresponding JSON coordinate files"

setup(
    name="LCT-utils",
    version=version,
    description=description,
    install_requires=[
        "numpy",
        "tifffile",
        'tqdm'
    ],
    author="Richard Qiu",
    packages=["LCT_utils"],
    entry_points={
        "console_scripts": [
            "crop-tiff=LCT_utils.crop_tiffs:main",
            "crop-json=LCT_utils.crop_jsons:main"
        ]}
)