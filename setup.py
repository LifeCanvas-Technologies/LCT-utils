from setuptools import setup

version = "0.2.1"
description = "Command-line cropping of TIFF files and stacks as well as corresponding JSON coordinate files"

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
        "console_scripts": [
            "crop-tiff=crop_tiffs.crop_tiffs:main",
            "crop-json=crop_tiffs.crop_jsons:main"
        ]}
)
