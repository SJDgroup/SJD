# SJD: R package for Structured Joint Decomposition

<!--
[![Travis build status](https://travis-ci.com/CHuanSite/SJD.svg?branch=main)](https://travis-ci.com/CHuanSite/SJD)
[![R-CMD-check](https://github.com/CHuanSite/SJD/workflows/R-CMD-check/badge.svg)](https://github.com/CHuanSite/SJD/actions)
-->

## Overview:
Structured Joint Decomposiiton (SJD) is a robust R package based on low rank models that provides novel decomposition methods for gene expression datasets that share similar biological progresses. It allows integrated dimensionality reduction for multiple datasets regardless of species, such as human and mouse datasets studying brain neuronal development.

SJD provides four decomposition algorithms: 1) separate decomposition; 2) concatenate decomposition; 3) joint decomposition; 4) two-stage sequential decomposition.

The first three algorithms allow three options of statistical computing process: Principal Component Analysis (PCA), Independent Component Analysis (ICA) and Nonnegative Matrix Factorization (NMF).
The last algorithm two-stage sequential decomposition is solely based on Principal Component Analysis (PCA).

All SJD algorithms require the same input data format: dataset_list, group and comp_num

- dataset_list: expression matrices to be analyzed
- group: the structure information of datasets to be analyzed
- comp_num: desired number of dimensionality components after decomposition
- weighting(optional): parameter specifying weights of dataset(s)

To install this package in R, run the following commands:

```R
library(devtools)
install_github("CHuanSite/SJD")
```

## Example Use

Please visit the vignette page [here](https://chuansite.github.io/SJD/)

[Structral Joint Decomposition (SJD)](https://chuansite.github.io/SJD/articles/StructralJointDecomposition.html)

[Two Stage Linked Component Analysis (2s-LCA)](https://chuansite.github.io/SJD/articles/twoStageLCA.html)

## Use SJD in Python

SJD can be used in Python by using a Python package called rpy2. Once rpy2 is installed in Python, [this](pythonExamples/SJD_in_python_example.ipynb) example can be followed as a guide.  
## Source

The latest developer version is available in the [Github repository](https://github.com/CHuanSite/SJD).

## Contact

For any improvements and issues that need to be addressed, please contact [Huan Chen](mailto:hzchenhuan@gmail.com) or [Shreyash Sonthalia](mailto:shreyash.gradschool@gmail.com).
