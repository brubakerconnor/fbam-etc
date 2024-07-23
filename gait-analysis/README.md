Analysis of Stride Interval Series
================

This directory contains the code to perform the analysis of stride
interval series collected by Hausdorff et al. (1998). To begin, run the
`gaitdata_download_process_fit.R` script after ensuring that the working
directory is the parent directory of the `gait-analysis` sub-directory.

``` sh
cd /path/to/fbam-etc/
Rscript gait-analysis/gaitdata_download_process_fit.R
```

This script will download the data directly from
[physionet.org](https://physionet.org/content/gaitndd/1.0.0/) and
process it according to the procedure detailed in the article. The
original data will be placed in the `gait-analysis/data-raw/` directory
and the output of this script will be placed in the
`gait-analysis/data/` directory.

## References

Hausdorff, J. M., Cudkowicz, M. E., Firtion, R., Wei, J. Y., &
Goldberger, A. L. (1998). *Gait variability and basal ganglia disorders:
Stride‐to‐stride variations of gait cycle timing in parkinson’s disease
and Huntington’s disease.* Movement Disorders, 13(3), 428–437.
<https://doi.org/10.1002/mds.870130310>
