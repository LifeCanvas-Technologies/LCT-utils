"""A module to generate cell count heat maps on annotated atlases
"""
import csv
import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tifffile
from skimage.segmentation import mark_boundaries, find_boundaries


def argparse(args=sys.argv[1:]):


def main(args=sys.argv[1:]):


def make_heatmap_fig(annotation,
                     cell_counts,
                     zs,
                     fig_layout,
                     figsize=[8,8],
                     cmap="hot_r",
                     round_to=10,
                     include_box=False):
