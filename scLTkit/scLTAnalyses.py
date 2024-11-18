from .utils import *

class LTAnalyses:
    def __init__(self, data_pre, data_pos, datasetName, sampleName, dataPath, savePath, lineage_identity,
                 pre_name="pre", pos_name="post", sel_cluster_name="cluster", cls_resolutions=None):
        print("------0. Preparing Basic information------")

        self.data_pre = data_pre
        self.data_pos = data_pos
        self.datasetName = datasetName
        self.sampleName = sampleName
        self.dataPath = dataPath
        self.savePath = savePath
        self.run_label = self.datasetName + "_" + self.sampleName
        self.run_label_time = self.run_label + '-' + pre_name + '_' + pos_name
        self.lineage_identity = lineage_identity
        self.pre_name = pre_name
        self.pos_name = pos_name
        self.sel_cluster_name = sel_cluster_name

        self.pre_colors = ["#43D9FE", "#E78AC3", "#FEC643", "#A6D854", "#FE6943", "#E5C494", "#33AEB1", "#FFEC1A",
                           "#4878D0", '#984EA3', '#CF747A', '#4DAF4A', '#C2C2C2', '#B5A2E0', '#F9B475']
        self.pos_colors = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#4878D0",
                           "#6ACC64", "#82C6E2", "#D65F5F", "#43D9FE", '#B383C1', '#BC8967', "#B3B3B3"]

        if cls_resolutions is None:
            cls_resolutions = [0.2, 0.2]
        self.cls_resolutions = cls_resolutions
        self.special_case = "Missing"

        self.barcodes_pre = None
        self.barcodes_pos = None
        self.cross_lin_mat = None
        self.diversity_metric = None
        self.flow_info = None
        self.all_de_df = None

        self.barcodes_pre, self.barcodes_pos = list(self.data_pre.obs[self.lineage_identity]), list(self.data_pos.obs[self.lineage_identity])
        self.cross_lin_mat = getLineageMatrix(bars=self.barcodes_pre, bars2=self.barcodes_pos)

        self.fig_runClonalHeterogeneity = []
        self.fig_runCellDynamics = None
        self.fig_runCelFateDiversity = None

        print("------End of prepareBasicInfo------")


    def runClonalHeterogeneity(self, method='Pearson', mode="cross"):
        print("------Mode: "+mode+" time-point------")
        if mode == "cross":
            sim_mat = getSimilarityMatrix(self.data_pre, self.data_pos, method=method)
            lin_mat = self.cross_lin_mat
            title = self.pre_name + ' -> ' + self.pos_name
        elif mode == "pre":
            sim_mat = getSimilarityMatrix(self.data_pre, self.data_pre, method=method)
            lin_mat = getLineageMatrix(bars=self.barcodes_pre, bars2=self.barcodes_pre)
            title = self.pre_name
        elif mode == "pos":
            sim_mat = getSimilarityMatrix(self.data_pos, self.data_pos, method=method)
            lin_mat = getLineageMatrix(bars=self.barcodes_pos, bars2=self.barcodes_pos)
            title = self.pos_name
        else:
            print("------Invalid mode selection for runClonalHeterogeneity------")
            return

        fig = plotSimilarityCompare(cross_sim=sim_mat,
                                    cross_lin_mat=lin_mat,
                                    title=title,
                                    savePath=self.savePath + self.run_label_time + '-SimDistrComp-hist-' + mode + '.png')

        self.fig_runClonalHeterogeneity.append(fig)
        del(sim_mat, lin_mat)
        gc.collect()
        print("------End of Mode: " + mode + " time-point------")


    def runCellDynamics(self):
        scd_obj = scStateDynamics(data_pre=self.data_pre, data_pos=self.data_pos,
                                  pre_name=self.pre_name, pos_name=self.pos_name,
                                  cls_prefixes=['', ''], run_label=self.run_label_time,
                                  pre_colors=self.pre_colors, pos_colors=self.pos_colors,
                                  savePath=self.savePath, saveFigFormat="png")
        scd_obj = runClustering(scd_obj, cls_resolutions=self.cls_resolutions)

        pre_fractions, pos_fractions = calculateFractions(scd_obj)
        t_row_sum = self.cross_lin_mat.sum(axis=1, keepdims=True)
        t_row_sum[t_row_sum == 0] = 1
        cls_lineage_mat = self.cross_lin_mat / t_row_sum
        cls_lineage_mat = np.array([[np.sum(cls_lineage_mat[np.where(scd_obj.data_pre.obs['cluster'] == str(i))[0]][:,
                                            np.where(scd_obj.data_pos.obs['cluster'] == str(j))[0]])
                                     for j in range(scd_obj.n_clus[1])] for i in range(scd_obj.n_clus[0])])
        cls_lineage_mat = cls_lineage_mat / scd_obj.data_pre.shape[0]
        flow_info = {'s': [str(i) for i in range(scd_obj.n_clus[0]) for j in range(scd_obj.n_clus[1])],
                     't': [str(j) for i in range(scd_obj.n_clus[0]) for j in range(scd_obj.n_clus[1])],
                     's_pm': list(cls_lineage_mat.reshape((1, -1))[0]),
                     't_pm': list(cls_lineage_mat.reshape((1, -1))[0])}
        flow_info = pd.DataFrame(flow_info)
        fig = plotFlowSankey(flow_info, self.pre_colors, self.pos_colors, pre_fractions=pre_fractions,
                             pos_fractions=pos_fractions,
                             figwidth=3.62, figheight=6, label_size=18, title=self.pre_name + '->' + self.pos_name,
                             label_position='twoside')
        fig.savefig(self.savePath + self.run_label_time + "-FlowSankey-lineage.png", dpi=300)
        self.fig_runCellDynamics = fig

        self.data_pre, self.data_pos, self.flow_info = scd_obj.data_pre, scd_obj.data_pos, flow_info


    def plotCellDynamicUMAP(self, fate_colname="Lineage_fate_label"):
        fig_size = 4
        pre_celltypes = np.unique(self.data_pre.obs[self.sel_cluster_name])
        fig, axs = plt.subplots(2+len(pre_celltypes), 1, figsize=(fig_size, fig_size*(2+len(pre_celltypes))))

        with plt.rc_context({'figure.figsize': (fig_size, fig_size)}):
            sc.pl.umap(self.data_pre, color=self.sel_cluster_name, palette=self.pre_colors, show=False, ax=axs[0])
            axs[0].set_title(self.pre_name)
            axs[0].axis('off')

        for i in range(len(pre_celltypes)):
            pre_celltype = pre_celltypes[i]
            temp_expr = self.data_pre.copy()
            temp_expr.obs[fate_colname] = temp_expr.obs[fate_colname].tolist()
            temp_expr.obs[fate_colname][self.data_pre.obs[self.sel_cluster_name] != pre_celltype] = "OtherClusters"
            temp_expr.obs[fate_colname] = temp_expr.obs[fate_colname].astype('category')

            original_categories = temp_expr.obs[fate_colname].unique().tolist()
            filtered_categories = [cat for cat in original_categories if cat not in [self.special_case, 'OtherClusters']]
            desired_order = filtered_categories + [self.special_case, 'OtherClusters']
            temp_expr.obs[fate_colname] = pd.Categorical(temp_expr.obs[fate_colname], categories=desired_order, ordered=True)

            pos_fates = np.array([s.split(' -> ')[1] for s in np.unique(temp_expr.obs[fate_colname]) if '->' in s])
            pos_celltype_uni = np.unique(self.data_pos.obs[self.sel_cluster_name])
            color_idx = [np.where(pos_celltype_uni == item)[0][0] for item in pos_fates]

            # color_idx = [int(s[-1]) for s in np.unique(temp_expr.obs[fate_colname]) if str(s[-1]).isdigit()]
            temp_color = [self.pos_colors[i] for i in color_idx]
            with plt.rc_context({'figure.figsize': (fig_size, fig_size)}):
                fate_colors = getColorMap(temp_expr.obs[fate_colname],
                                          [self.special_case, 'OtherClusters'],
                                          default_palette=temp_color)
                sc.pl.umap(temp_expr, color=fate_colname, palette=fate_colors, show=False, ax=axs[1+i])
                axs[1+i].set_title('Cell Fate of '+self.pre_name+'-cluster' + str(pre_celltype))
                axs[1+i].axis('off')

        with plt.rc_context({'figure.figsize': (fig_size, fig_size)}):
            sc.pl.umap(self.data_pos, color=self.sel_cluster_name, palette=self.pos_colors, show=False, ax=axs[1+len(pre_celltypes)])
            axs[1+len(pre_celltypes)].set_title(self.pos_name)
            axs[1+len(pre_celltypes)].axis('off')

        # plt.tight_layout()
        plt.savefig(self.savePath + self.run_label_time + '_combined_umap.png', dpi=300, bbox_inches='tight')
        plt.show()
        self.fig_runCelFateDiversity = fig


    def runCellFateDiversity(self):

        pos_celltypes = self.data_pos.obs[self.sel_cluster_name]
        cell_2lin_cls = np.array([(self.cross_lin_mat[:, np.where(pos_celltypes == celltype)[0]] > 0).sum(axis=1).tolist()
                                  for celltype in np.unique(pos_celltypes)]).T

        fate_cls = np.unique(pos_celltypes)[np.argmax(cell_2lin_cls, axis=1)]
        fate_cls[np.sum(cell_2lin_cls, axis=1) == 0] = self.special_case
        self.data_pre.obs['Lineage_fate'] = fate_cls

        fate_cls = self.data_pre.obs[self.sel_cluster_name].astype(str) + ' -> ' + fate_cls
        fate_cls[np.sum(cell_2lin_cls, axis=1) == 0] = self.special_case
        self.data_pre.obs['Lineage_fate_label'] = fate_cls

        pre_df = self.data_pre.obs[[self.sel_cluster_name, 'Lineage_fate', 'Lineage_fate_label']]
        pre_df.to_csv(self.savePath + self.run_label_time + '-lineage-analysis' + '.txt', sep='\t')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.plotCellDynamicUMAP()

        cfrs, afs = fateVectorDiversity(self.data_pre, cell_2lin_cls, fate_cls_name="Lineage_fate")
        ncs, ecs = neighborFateLabelDiversity(self.data_pre, fate_cls_name="Lineage_fate")
        self.diversity_metric = [ncs, ecs, cfrs, afs]

        print("Neighboring cell fate label consistency (NFC): {:.4f}".format(ncs))
        print("Neighboring cell fate similarity (NFS): {:.4f}".format(afs))
        print("Cell fate randomness (CFR): {:.4f}".format(cfrs))
        print("Neighboring cell fate label randomness (NFR) : {:.4f}".format(ecs))


    def runSubClusterDiff(self, sel_cls=None, sel_fates=None):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if sel_fates is None:
                all_de_df = pd.DataFrame()
                print("Mode of runSubClusterDiff: 1 v.s. Rest")
                fate_colname = 'Lineage_fate'
                savename = "_DE_all_fate_genes"
                pre_celltypes = np.unique(self.data_pre.obs[self.sel_cluster_name])
                # Select one cluster in pre-timepoint, use 'Lineage_fate' for diff analysis
                for si in pre_celltypes:
                    cur_expr = self.data_pre[self.data_pre.obs[self.sel_cluster_name] == str(si)]
                    cur_expr = cur_expr[cur_expr.obs[fate_colname] != self.special_case]
                    t_v, t_n = np.unique(cur_expr.obs[fate_colname], return_counts=True)
                    if np.all(t_n > 1) == False:
                        mask = np.isin(np.array(cur_expr.obs[fate_colname]), t_v[np.argwhere(t_n > 1).flatten()])
                        cur_expr = cur_expr[np.where(mask)[0]]
                        t_v, t_n = np.unique(cur_expr.obs[fate_colname], return_counts=True)
                    try:
                        sc.tl.rank_genes_groups(cur_expr, groupby=fate_colname, method='wilcoxon')
                        # Select genes
                        for i, f in enumerate(cur_expr.obs[fate_colname].unique()):
                            fate_str = str(si) + '->' + str(f)
                            de_df = generateDEGs(cur_expr, index=i, cluster=si, fate_str=fate_str)
                            all_de_df = pd.concat([all_de_df, de_df], axis=0)

                    except ZeroDivisionError as e:
                        pass

                self.all_de_df = all_de_df

                if all_de_df.shape[0] > 0:
                    print(all_de_df)
                    all_de_df.to_csv(self.savePath + self.run_label_time + savename + '-onlyLT.txt',
                                     sep='\t', header=True, index=False)

            else:
                print("Mode of runSubClusterDiff: 1 v.s. 1")
                fate_colname = 'Lineage_fate_label'

                for idx in range(len(sel_fates)):
                    sel_fate = sel_fates[idx]
                    print(sel_fate)
                    savename = '_DE_fate_genes-' + sel_cls + str(idx+1)
                    all_de_df = pd.DataFrame()
                    # Only one cluster at pre-timepoint, use 'Lineage_fate_label' for diff analysis
                    cur_expr = self.data_pre[self.data_pre.obs[fate_colname].isin(sel_fate)]
                    t_v, t_n = np.unique(cur_expr.obs[fate_colname], return_counts=True)
                    if np.all(t_n > 1) == False:
                        mask = np.isin(np.array(cur_expr.obs[sel_fate]), t_v[np.argwhere(t_n > 1).flatten()])
                        cur_expr = cur_expr[np.where(mask)[0]]
                        t_v, t_n = np.unique(cur_expr.obs[sel_fate], return_counts=True)
                    try:
                        # cur_expr.uns['log1p']["base"] = None
                        sc.tl.rank_genes_groups(cur_expr, groupby=fate_colname, method='wilcoxon')
                        subfates = [name for name, _ in cur_expr.uns['rank_genes_groups']['names'].dtype.fields.items()]
                        for i in range(len(subfates)):
                            de_df = generateDEGs(cur_expr, index=i, cluster=sel_cls, fate_str=subfates[i])
                            all_de_df = pd.concat([all_de_df, de_df], axis=0)

                    except ZeroDivisionError as e:
                        pass

                    if all_de_df.shape[0] > 0:
                        all_de_df.to_csv(self.savePath + self.run_label_time + savename + '-onlyLT.txt',
                                         sep='\t', header=True, index=False)

                    print(all_de_df)


    def runLTAnalyses(self):
        print("------1. Start of runClonalHeterogeneity------")
        self.runClonalHeterogeneity(method='Pearson', mode="cross")
        self.runClonalHeterogeneity(method='Pearson', mode="pre")
        self.runClonalHeterogeneity(method='Pearson', mode="pos")
        print("------End of runClonalHeterogeneity------")

        print("------2. Start of runCellDynamics------")
        self.runCellDynamics()
        print("------End of runCellDynamics------")

        print("------3. Start of runCellFateDiversity------")
        self.runCellFateDiversity()
        print("------End of runCellFateDiversity------")

        print("------4. Start of runSubClusterDiff------")
        self.runSubClusterDiff()
        print("------End of runSubClusterDiff------")
