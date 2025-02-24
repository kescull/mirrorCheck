# mirrorCheck
## Facilitator functions for getting and assessing DESeq2 lfcShrink results
This package has just two functions. The first, `run_DESeq_all_contrasts()`, is designed to help [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) users extract results for all the possible pairwise comparisons after running `DESeq()` on a DESeqDataSet with multiple groups to compare, as a wrapper for calling either `results()` or `lfcShrink()` iteratively. If `lfcShrink()` was chosen, the second function, `compare_reciprocal_contrasts()` can help to visualise how well 
`lfcShrink()` was able to perform with the given data, via diagnostic plots which analyse how well the pairs of reciprocal contrasts agree (i.e. when using alternative reference groups).

### Installation
You can install this package from GitHub using the following code in R (first installing `devtools` if necessary):
```
# install.packages("devtools")
devtools::install_github("kescull/mirrorCheck")
```
### Citation
If used, please cite our paper (...)


