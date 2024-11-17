# scLT-kit

## Introduction
`scLTkit` is a toolkit for analyzing single-cell lineage-tracing (LT-scSeq) data.

## System Requirements
- Python version: >= 3.7

## Dependencies
- numpy, pandas, sklearn, scipy, matplotlib, seaborn

## Installation

The Release version `scLTkit` python package will be published soon. It can be installed directly via pip:
```
pip install scLTkit
```

## Quick Start of LT-scSeq data enhancement

Refer to folder: [tutorial](https://github.com/czythu/scTrace/tree/main/tutorial) for full pipeline.

Example data: [Larry-Invitro-differentiation](https://cloud.tsinghua.edu.cn/f/1b94b3229f4a4c52985e/?dl=1)

Below are the introduction to the main steps in `scLTkit`.

1. `runLTStatistics`: Main steps of Class `LTStatistics`, including `getBarcodingFractions` and `getClonalSizes`.

2. `runLTAnalyses`: Main steps of Class `LTAnalyses`, including `runClonalHeterogeneity` (within & cross time-point),
`runCellDynamics`, `runCellFateDiversity`, and `runSubClusterDiff`


