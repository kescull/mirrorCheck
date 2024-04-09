# lfcShrinkCheck
## Facilitator functions for getting and assessing DESeq2 lfcShrink results
This package has just two functions. The first, `run.DESeq.all.contrasts()`, is designed to help DESeq2 users extract all the possible pairwise comparisons after running `DESeq()` on a DESeqDataSet possessing a simple design, 
such as `design = ~condition`, as a wrapper for calling either `results()` or `lfcShrink()` (using apeglm algorithm) iteratively. If `lfcShrink()` was chosen, the second function, `compare.reciprocal.contrasts()` can help to visualise how well 
`lfcShrink()` was able to perform with the given data, via diagnostic plots which analyse how well the pairs of reciprocal contrasts agree (i.e. when using alternative reference groups).

### Installation
You can install this package from GitHub using the following code in R (first installing `devtools` if necessary):
```
# install.packages("devtools")
devtools::install_github("kescull/lfcShrinkCheck")
```
### Citation
If used, please cite our paper (...)


