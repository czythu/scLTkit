from .utils import *

class LTStatistics:
    def __init__(self, data_pre, data_pos, datasetName, sampleName, dataPath, savePath, lineage_identity):
        print("------0. Preparing Basic information------")

        self.data_pre = data_pre
        self.data_pos = data_pos
        self.datasetName = datasetName
        self.sampleName = sampleName
        self.dataPath = dataPath
        self.savePath = savePath
        self.run_label = self.datasetName + "_" + self.sampleName
        self.lineage_identity = lineage_identity

        self.barcodes_pre = None
        self.barcodes_pos = None
        self.pre_barcode_den = None
        self.pos_barcode_den = None
        self.flow_out_den = None
        self.flow_in_den = None

        self.size_freq_pre = None
        self.size_freq_pos = None

        self.fig_getBarcodeFractions = None
        self.fig_getClonalSizes = None
        print("------End of prepareBasicInfo------")


    def getBarcodingFractions(self):
        self.barcodes_pre, self.barcodes_pos = list(self.data_pre.obs[self.lineage_identity]), list(self.data_pos.obs[self.lineage_identity])
        # Within-timepoint & Cross-timepoint
        cross_lin_mat = getLineageMatrix(bars=self.barcodes_pre, bars2=self.barcodes_pos)
        n_pre, n_pre_nan, n_out = (cross_lin_mat.shape[0], Counter(self.barcodes_pre)[np.nan], sum(cross_lin_mat.sum(axis=1) != 0))
        n_pos, n_pos_nan, n_in = (cross_lin_mat.shape[1], Counter(self.barcodes_pos)[np.nan], sum(cross_lin_mat.sum(axis=0) != 0))
        self.pre_barcode_den, self.pos_barcode_den = 1 - n_pre_nan/n_pre, 1 - n_pos_nan/n_pos
        self.flow_out_den, self.flow_in_den = n_out/n_pre, n_in/n_pos

        print("------Pre time point------")
        print("Number of cells in the former time point: ", n_pre)
        print("Number of cells with lineage barcode: ", n_pre-n_pre_nan)
        print("Number of cells with flow-out information: ", n_out)
        print("Barcoding fraction of pre-timepoint: {:.4f}".format(self.pre_barcode_den))
        print("Flow-out density of pre-timepoint: {:.4f}".format(self.flow_out_den))
        print("------Pos time point------")
        print("Number of cells in the latter time point: ", n_pos)
        print("Number of cells with lineage barcode: ", n_pos-n_pos_nan)
        print("Number of cells with flow-in information: ", n_in)
        print("Barcoding fraction of pos-timepoint: {:.4f}".format(self.pos_barcode_den))
        print("Flow-in density of pos-timepoint: {:.4f}".format(self.flow_in_den))

        self.fig_getBarcodeFractions = plotBarcodeFraction([n_pre_nan, n_pre-n_pre_nan-n_out, n_out],
                                                           [n_pos_nan, n_pos-n_pos_nan-n_in, n_in],
                                                           savePath=self.savePath + self.run_label + '-BarcodeFrac.png')


    def getClonalSizes(self):
        clone_size_pre = pd.Series(self.barcodes_pre).dropna().value_counts()
        clone_size_pos = pd.Series(self.barcodes_pos).dropna().value_counts()
        self.size_freq_pre = clone_size_pre.value_counts()
        self.size_freq_pos = clone_size_pos.value_counts()

        self.fig_getClonalSizes = plotClonalSizes(size_freq_pre=self.size_freq_pre,
                                                  size_freq_pos=self.size_freq_pos,
                                                  savePath=self.savePath + self.run_label + '-ClonalSizes.png')


    def runLTStatistics(self):
        print("------1. Start of getBarcodingFractions------")
        self.getBarcodingFractions()
        print("------End of getBarcodingFractions------")

        print("------2. Start of getClonalSizes------")
        self.getClonalSizes()
        print("------End of getClonalSizes------")
