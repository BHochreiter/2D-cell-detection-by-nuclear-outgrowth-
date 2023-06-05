# 2D-cell-detection-by-nuclear-outgrowth

Created by Bernhard Hochreiter 2021

Published unter GNU General Public License 3.0

Cells in close proximity are often impossible to segment for cell-based measurements. This imageJ macro can approximate individual cells via an outgrowth algorithm from the nucleus.

First, a nuclear stain is detected by intensity thresholding and the single nuclei are detected as objects by size exclusion. In order to measure cytosolic intensity, an outgrowth algorithm grows each detected nuclei by adding the pixels in direct proximity. When two cells would connect during this step, the algorithm stops, keeping single cells separated. This is repeated multiple times, with the number of successive outgrowth iterations set by the user.

The output consists of 1) an image displaying the detected and outgrown area of each individual cell either in monochrome or rainbow color and 2) a list of average intensities in nucleus, cytosol and total for each cell and all supplied channels in the dataset.
