# Deep-learning-based-intracellular-foci-enrichment-analysis
# The deep learning image segmentation algorithm cellpose (Python) is used to produce detailed cell segmentation. 
# The segmented images are used in turn to identify (MATLAB) foci in a target channel and calculate the intensity of signals surrounding the foci in other channels. The output is plotted as signal intensity over distance from the foci perimeter pixels. 
# Cellpose segmentation was implemented using the command line option as described in https://cellpose.readthedocs.io/en/latest/command.html
# Cellpose was developed by Carsen Stringer and Marius Pachitariu and presented in the seminal paper "Cellpose: a generalist algorithm for cellular segmentation, Nature Methods doi: 10.1038/s41592-020-01018-x'
