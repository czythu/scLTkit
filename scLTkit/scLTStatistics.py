from .utils import *

class LTStatistics:
    """
    A class for performing lineage tracing statistics, including barcoding fractions
    and clonal size distributions.

    Attributes:
        data_pre: Data from the pre-timepoint.
        data_pos: Data from the post-timepoint.
        datasetName: Name of the dataset being analyzed.
        sampleName: Name of the sample being analyzed.
        dataPath: Path to input data.
        savePath: Path to save outputs.
        lineage_identity: Column name representing lineage identity.
        barcodes_pre: List of barcodes from the pre-timepoint.
        barcodes_pos: List of barcodes from the post-timepoint.
        pre_barcode_den: Barcoding fraction at the pre-timepoint.
        pos_barcode_den: Barcoding fraction at the post-timepoint.
        flow_out_den: Flow-out density from the pre-timepoint.
        flow_in_den: Flow-in density at the post-timepoint.
        size_freq_pre: Frequency distribution of clone sizes at the pre-timepoint.
        size_freq_pos: Frequency distribution of clone sizes at the post-timepoint.
        fig_getBarcodeFractions: Figure of barcoding fractions.
        fig_getClonalSizes: Figure of clonal size distributions.
    """

    def __init__(self, data_pre, data_pos, datasetName, sampleName, dataPath, savePath, lineage_identity):
        """
        Initializes the LTStatistics class with basic information and dataset attributes.

        Args:
            data_pre: Pre-timepoint data.
            data_pos: Post-timepoint data.
            datasetName: Name of the dataset.
            sampleName: Name of the sample.
            dataPath: Path to the dataset.
            savePath: Path to save results.
            lineage_identity: Column name for lineage identity in the data.
        """
        print("------0. Preparing Basic information------")

        self.data_pre = data_pre
        self.data_pos = data_pos
        self.datasetName = datasetName
        self.sampleName = sampleName
        self.dataPath = dataPath
        self.savePath = savePath
        self.run_label = self.datasetName + "_" + self.sampleName
        self.lineage_identity = lineage_identity

        # Initialize attributes for lineage statistics
        self.barcodes_pre = None
        self.barcodes_pos = None
        self.pre_barcode_den = None
        self.pos_barcode_den = None
        self.flow_out_den = None
        self.flow_in_den = None

        self.size_freq_pre = None
        self.size_freq_pos = None

        # Initialize figures
        self.fig_getBarcodeFractions = None
        self.fig_getClonalSizes = None
        print("------End of prepareBasicInfo------")

    def getBarcodingFractions(self):
        """
        Calculates the barcoding fractions and flow densities at the pre- and post-timepoints.
        Generates a visualization of the barcoding fractions.
        """
        # Extract barcodes from lineage identity column
        self.barcodes_pre, self.barcodes_pos = list(self.data_pre.obs[self.lineage_identity]), list(self.data_pos.obs[self.lineage_identity])

        # Compute cross-lineage matrix between pre- and post-timepoints
        cross_lin_mat = getLineageMatrix(bars=self.barcodes_pre, bars2=self.barcodes_pos)

        # Calculate barcoding and flow densities
        n_pre, n_pre_nan, n_out = (cross_lin_mat.shape[0], Counter(self.barcodes_pre)[np.nan], sum(cross_lin_mat.sum(axis=1) != 0))
        n_pos, n_pos_nan, n_in = (cross_lin_mat.shape[1], Counter(self.barcodes_pos)[np.nan], sum(cross_lin_mat.sum(axis=0) != 0))
        self.pre_barcode_den, self.pos_barcode_den = 1 - n_pre_nan / n_pre, 1 - n_pos_nan / n_pos
        self.flow_out_den, self.flow_in_den = n_out / n_pre, n_in / n_pos

        # Print results for debugging and reference
        print("------Pre time point------")
        print("Number of cells in the former time point: ", n_pre)
        print("Number of cells with lineage barcode: ", n_pre - n_pre_nan)
        print("Number of cells with flow-out information: ", n_out)
        print("Barcoding fraction of pre-timepoint: {:.4f}".format(self.pre_barcode_den))
        print("Flow-out density of pre-timepoint: {:.4f}".format(self.flow_out_den))
        print("------Pos time point------")
        print("Number of cells in the latter time point: ", n_pos)
        print("Number of cells with lineage barcode: ", n_pos - n_pos_nan)
        print("Number of cells with flow-in information: ", n_in)
        print("Barcoding fraction of pos-timepoint: {:.4f}".format(self.pos_barcode_den))
        print("Flow-in density of pos-timepoint: {:.4f}".format(self.flow_in_den))

        # Generate and save visualization of barcoding fractions
        self.fig_getBarcodeFractions = plotBarcodeFraction([n_pre_nan, n_pre - n_pre_nan - n_out, n_out],
                                                           [n_pos_nan, n_pos - n_pos_nan - n_in, n_in],
                                                           savePath=self.savePath + self.run_label + '-BarcodeFrac.png')

    def getClonalSizes(self):
        """
        Calculates the clonal size distributions at the pre- and post-timepoints.
        Generates a visualization of the clonal size distributions.
        """
        # Compute clone sizes by counting occurrences of barcodes
        clone_size_pre = pd.Series(self.barcodes_pre).dropna().value_counts()
        clone_size_pos = pd.Series(self.barcodes_pos).dropna().value_counts()

        # Compute frequency distributions of clone sizes
        self.size_freq_pre = clone_size_pre.value_counts()
        self.size_freq_pos = clone_size_pos.value_counts()

        # Generate and save visualization of clonal size distributions
        self.fig_getClonalSizes = plotClonalSizes(size_freq_pre=self.size_freq_pre,
                                                  size_freq_pos=self.size_freq_pos,
                                                  savePath=self.savePath + self.run_label + '-ClonalSizes.png')

    def runLTStatistics(self):
        """
        Executes the full lineage tracing statistics pipeline:
        1. Calculates barcoding fractions.
        2. Analyzes clonal size distributions.
        """
        print("------1. Start of getBarcodingFractions------")
        self.getBarcodingFractions()
        print("------End of getBarcodingFractions------")

        print("------2. Start of getClonalSizes------")
        self.getClonalSizes()
        print("------End of getClonalSizes------")
