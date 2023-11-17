# lsd_fmri_fc_entropy

1. Applies connReader on a CONN project to extract network-level functional connectivity values.
    input: ROI CONN project, output: csv file in long format with functional connectivity values at network-level
2. Computes entropy and graph theory metrics (including normalized modularity) using CopBET and BCT toolboxes
    - for BOLD complexity, DCC and meta-state complexity, input: ROI denoised timeseries
    - for path-length distribution entropy and graph theory metrics, input: ROI-ROI connectivity matrix
    - for multi-scale sample entropy, input: voxel denoised timeseries
3. Add entropy and graph theory metrics to a csv file in long format
4. Compute statistics with linear mixed effect models and plot fonctional connectivity matrices, normalized modularity and entropy metrics
    input: csv file, output: plots saved as png and svg, and tables with statistics
5. Transform the data csv from a long format to a wide format
    input: csv files as long format, output: data_wide.csv
