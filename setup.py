from setuptools import setup

version = "0.3.0"
description = "Command-line cropping of TIFF files and stacks as well as corresponding JSON coordinate files"

setup(
    name="LCT-utils",
    version=version,
    description=description,
    install_requires=[
        "numpy",
        "scipy",
        "tifffile",
        'tqdm',
        'scikit-image'
    ],
    author="Richard Qiu",
    packages=["LCT_utils"],
    entry_points={
        "console_scripts": [
            "transform-tiff=LCT_utils.transform_tiffs:main",
            "transform-json=LCT_utils.transform_jsons:main",
            "json-to-tiff=LCT_utils.json_to_tiff:main"
        ]}
)
