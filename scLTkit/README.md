## scLTkit function Descriptions
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