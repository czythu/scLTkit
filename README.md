# scLT-kit

[![PyPI version](https://badge.fury.io/py/scLTkit.svg?icon=si%3Apython)](https://pypi.org/project/scLTkit/)

## Introduction
`scLT-kit` is a toolkit for analyzing single-cell lineage-tracing (LT-scSeq) data.

<div align=center>
<img alt="scLT-kit workflow" src="tutorial/scLTkit.png" width="90%"/>
</div>

## System Requirements
- Python version: >= 3.7

## Dependencies
- numpy, pandas, sklearn, scipy, matplotlib, seaborn

## Installation

The release version of `scLT-kit`+ python package can be installed directly via pip:
```
pip install scLTkit
```

## Quick Start of scLT-kit & Example datasets

Refer to folder: [tutorial](https://github.com/czythu/scLTkit/tree/main/tutorial) for full pipeline.

Example data1: [Larry-Invitro-differentiation](https://cloud.tsinghua.edu.cn/f/1b94b3229f4a4c52985e/?dl=1)

Example data2: [TraCe-seq-tumor](https://cloud.tsinghua.edu.cn/f/dae5b3ff8bd04177bd5f/?dl=1)


## File Descriptions
- scLTkit functions
  - [scLTkit/scLTStatistics.py](https://github.com/czythu/scLTkit/blob/main/scLTkit/scLTStatistics.py): This file contains the class `LTStatistics`,
  which includes three main functions:
    - `getBarcodingFractions` calculates the proportion of lineage barcodes, evaluating the coverage of the barcode labeling process in lineage tracing.
    - `getClonalSizes` computes the distribution of clonal sizes, visualizing the number of cells derived from each clone.
    - `runLTStatistics` is the "forward function" of class `LTStatistics`, including `getBarcodingFractions` and `getClonalSizes`.
  - [scLTkit/scLTAnalyses.py](https://github.com/czythu/scLTkit/blob/main/scLTkit/scLTAnalyses.py): This file contains the class `Analyses`,
  which provides five primary functions:
    - `runClonalHeterogeneity` analyzes the similarity between cells within a clone and across different clones, offering insights into clonal heterogeneity and lineage relationships.
    - `runCellDynamics` visualizes cluster-level cell fate flow through a Sankey plot.
    - `runCellFateDiversity` quantifies cell fate diversity using four evaluation metrics.
    - `runSubClusterDiff` performs differential expression analysis on subpopulations with different lineage trajectories, identifying key genes associated with distinct cell fates.
    - `runLTAnalyses` is the "forward function" of class `LTAnalyses`, including `runClonalHeterogeneity` (within & cross time-point), `runCellDynamics`, `runCellFateDiversity`, and `runSubClusterDiff`.
  - [scLTkit/utils.py](https://github.com/czythu/scLTkit/blob/main/scLTkit/utils.py): This file contains basic utility functions that support various aspects of data analysis,
  including clustering, metric calculation (e.g., transcriptomic similarity and fate diversity), differential analysis, and visualization.

- Tutorial files
  - [tutorial/Larry-InvitroDiff.ipynb](https://github.com/czythu/scLTkit/blob/main/tutorial/Larry-InvitroDiff.ipynb): This Jupyter notebook contains the complete analysis pipeline and results for the hematopoietic differentiation data
  presented in `Figure 1` of [Caleb Weinreb et al.](https://www.science.org/doi/10.1126/science.aaw3381) The analysis follows two key steps `runLTStatistics` and `runLTAnalyses`.
  Except for the differential analysis including all cell subtypes, the analysis are also performed specifically on several cell subtypes of interest.
  - [tutorial/TraCeseq-tumor.ipynb](https://github.com/czythu/scLTkit/blob/main/tutorial/TraCeseq-tumor.ipynb): This Jupyter notebook contains the complete analysis pipeline and results for the Erlotinib-perturbed lung cancer cell line
  presented in `Figure 1` of [Matthew T. Chang et al.](https://www.nature.com/articles/s41587-021-01005-3) The analysis also follows two key steps `runLTStatistics` and `runLTAnalyses`.


## Technical Details
Refer to folder: [tutorial](https://github.com/czythu/scLTkit/tree/main/tutorial/) for technical details.
