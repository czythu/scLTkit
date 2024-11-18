import warnings
import os
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import math
from sklearn.metrics import pairwise
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, rgb2hex
import seaborn as sns
from collections import Counter


# Modified from scStateDynamics
def runClustering(scd, cls_resolutions, savePath=None):
    print("## Identify cell clusters and cell states")
    if savePath is None:
        savePath = scd.savePath

    scd.cls_resolutions = cls_resolutions

    sc.tl.leiden(scd.data_pre, key_added='cluster', resolution=scd.cls_resolutions[0])
    sc.tl.leiden(scd.data_pos, key_added='cluster', resolution=scd.cls_resolutions[1])

    scd.data_pre.obs['cluster_id'] = scd.data_pre.obs['cluster']
    scd.data_pre.obs['cluster'] = [scd.cls_prefixes[0] + i for i in scd.data_pre.obs['cluster_id']]
    scd.data_pos.obs['cluster_id'] = scd.data_pos.obs['cluster']
    scd.data_pos.obs['cluster'] = [scd.cls_prefixes[1] + i for i in scd.data_pos.obs['cluster_id']]

    scd.n_clus = np.array([len(scd.data_pre.obs["cluster"].unique()), len(scd.data_pos.obs["cluster"].unique())])

    print('| - N_cells:', [scd.data_pre.shape[0], scd.data_pos.shape[0]])
    print('| - N_clusters:', scd.n_clus)

    saveName = savePath + scd.run_label + '_Cell-PairUMAP.' + scd.saveFigFormat

    scd.p_clusterUMAP = scd.plotPairScatter(value='cluster', axis_value='X_umap', x_lab='UMAP_1', y_lab='UMAP_2',
                                            saveName=saveName)

    return scd


# Modified from scStateDynamics
def plotScatter(sc_obj, value='cluster', axis_value='X_umap', labs=['UMAP_1', 'UMAP_2'],
                title=None, palette=None, legend_title='Cluster', saveFig=False, saveName='CellScatter.png'):
    p_data = pd.DataFrame({labs[0]: sc_obj.obsm[axis_value][:, 0],
                           labs[1]: sc_obj.obsm[axis_value][:, 1],
                           value: list(sc_obj.obs[value])})
    p_data = p_data.sort_values(by=value, axis=0, ascending=True)
    p_ratio = p_data[labs].max() - p_data[labs].min()
    p_ratio = p_ratio[labs[0]] / p_ratio[labs[1]]

    if saveFig:
        fig = plt.figure(figsize=(3.5, 3), tight_layout=True)
    ax = sns.scatterplot(data=p_data, x=labs[0], y=labs[1], hue=value, palette=palette,
                         s=5, alpha=0.6, linewidth=0)
    if legend_title is None:
        legend_title = value
    ax.legend(title=legend_title, markerscale=3, loc=6, bbox_to_anchor=(1.01, 0.5), ncol=1, handletextpad=0)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, loc='left', fontsize=13)
    ax.set_aspect(p_ratio)

    if saveFig:
        fig.savefig(saveName, dpi=300)
        print('| - Saving figure:', saveName)

    return (ax)


# Modified from scStateDynamics
class scStateDynamics:

    def __init__(self, data_pre, data_pos, pre_name, pos_name, run_label, pre_colors=None, pos_colors=None,
                 cls_prefixes=['S', 'T'], savePath="", saveFigFormat="png"):
        self.data_pre = data_pre
        self.data_pos = data_pos
        self.pre_name = pre_name
        self.pos_name = pos_name
        self.run_label = run_label
        self.savePath = savePath
        if not os.path.exists(savePath):
            os.makedirs(savePath)
        self.saveFigFormat = saveFigFormat

        self.pre_colors = pre_colors
        self.pos_colors = pos_colors
        if self.pre_colors is None:
            self.pre_colors = ["#FEC643", "#437BFE", "#43FE69", "#FE6943", "#E78AC3",
                               "#43D9FE", "#FFEC1A", "#E5C494", "#A6D854", "#33AEB1",
                               "#EA6F5D", "#FEE8C3", "#3DBA79", "#D6EBF9", "#7F6699",
                               "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                               "#333399", "#679966", "#c12e34", "#66d5a5", "#5599ff"]
        if self.pos_colors is None:
            self.pos_colors = ['#EE854A', '#4878D0', '#6ACC64', '#D65F5F', '#956CB4',
                               '#82C6E2', '#D5BB67', '#8C613C', '#DC7EC0', '#797979']
        self.cls_prefixes = cls_prefixes
        self.cls_resolutions = [0.5, 0.5]
        self.n_clus = np.array([])

        self.p_clusterUMAP = None
        self.bool_Clustering = False

    def __repr__(self):
        all_str = "scStateDynamics object: %s" % self.run_label

        if self.savePath != "":
            all_str = all_str + "\n+ 0. Data information: \n  - savePath: " + self.savePath

        all_str = all_str + "\n  - Pre-data dimensions: " + str(self.data_pre.shape)
        all_str = all_str + "\n  - Post-data dimensions: " + str(self.data_pos.shape)

        if self.bool_Clustering:
            all_str = all_str + "\n+ 1. Identifying cell clusters and cell states is finished."
            all_str = all_str + "\n  - Number of clusters: " + str(self.n_clus)

        return (all_str)

    ### Identify cell clusters and cell states
    def plotPairScatter(self, value, axis_value, x_lab, y_lab, pre_colors=None, pos_colors=None,
                        titles=None, saveFig=True, saveName='Cell-PairScatter.png'):
        if isinstance(value, str):
            value = [value, value]

        if pre_colors is None:
            pre_colors = self.pre_colors[0:len(self.data_pre.obs[value[0]].unique())]
        if pos_colors is None:
            pos_colors = self.pos_colors[0:len(self.data_pos.obs[value[1]].unique())]

        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(1, 2, 1)

        if titles is None:
            pre_title = self.pre_name
            pos_title = self.pos_name
        elif type(titles) == list:
            if len(titles) == 2:
                pre_title, pos_title = titles
        else:
            print('The parameter `titles` should be a list with two items. Use the default titles instead.')
            pre_title = self.pre_name
            pos_title = self.pos_name

        plotScatter(self.data_pre, value=value[0], axis_value=axis_value, labs=[x_lab, y_lab],
                    title=pre_title, palette=pre_colors)
        ax2 = fig.add_subplot(1, 2, 2)
        plotScatter(self.data_pos, value=value[1], axis_value=axis_value, labs=[x_lab, y_lab],
                    title=pos_title, palette=pos_colors)
        fig.tight_layout()

        if saveFig:
            fig.savefig(saveName, dpi=300)
            print('| - Saving figure:', saveName)

        return fig


def getSimilarityMatrix(data_pre, data_pos, method='Pearson'):
    data_comb = np.vstack((data_pre.raw.to_adata()[:, data_pre.var.index[data_pre.var.highly_variable]].X.toarray(),
                           data_pos.raw.to_adata()[:, data_pos.var.index[data_pos.var.highly_variable]].X.toarray()))
    n_pre = data_pre.shape[0]
    n_pos = data_pos.shape[0]

    all_sim = None
    if method == 'Pearson':
        all_sim = np.corrcoef(data_comb)
    elif method == 'Cosine':
        all_sim = 1 - pairwise.cosine_distances(data_comb)

    cross_sim = all_sim[0:n_pre, n_pre:(n_pre + n_pos)]
    return cross_sim


def getLineageMatrix(bars, bars2, ignore_bars=None):
    bars = pd.Series(bars, dtype='category')
    bars2 = pd.Series(bars2, dtype='category')
    lin_mat = pd.DataFrame(np.zeros((len(bars), len(bars2))))
    # Intersection of 2 timepoints
    comm_bars = list(set(bars).intersection(set(bars2)) - {np.nan})
    if ignore_bars:
        comm_bars = list(set(comm_bars) - set(ignore_bars))
    # Search in common barcodes
    for bar in comm_bars:
        t_ix1 = np.where(bars == bar)[0]
        t_ix2 = np.where(bars2 == bar)[0]
        lin_mat.iloc[t_ix1, t_ix2] = 1

    return np.array(lin_mat)


def plotBarcodeFraction(pre_info, pos_info, savePath):
    colors_pre = ["#B3E5FC", "#81D4FA", "#4FC3F7"]
    colors_pos = ["#E1C6E6", "#DDA9E4", "#A45DB8"] # "#BA68C8","#E4B3ED", EDE1F8
    #
    labels_pre = ['Missing', 'Only pre', 'Transfer to post']
    labels_pos = ['Missing', 'Only post', 'Inherited from pre']

    fig, axs = plt.subplots(1, 2, figsize=(9, 3))

    wedges1, texts1, autotexts1 = axs[0].pie(pre_info, colors=colors_pre, autopct='%1.1f%%', startangle=90)
    axs[0].axis('equal')
    axs[0].set_title('Pre-timepoint barcoding fraction')
    axs[0].legend(wedges1, labels_pre, title="Pre-barcoding", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    wedges2, texts2, autotexts2 = axs[1].pie(pos_info, colors=colors_pos, autopct='%1.1f%%', startangle=90)
    axs[1].axis('equal')
    axs[1].set_title('Post-timepoint barcoding fraction')
    axs[1].legend(wedges2, labels_pos, title="Post-barcoding", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    plt.tight_layout()
    plt.show()

    fig.savefig(savePath, dpi=300)

    return fig


def plotClonalSizes(size_freq_pre, size_freq_pos, savePath):
    fig, axs = plt.subplots(1, 2, figsize=(8, 3))

    sorted_counts1 = size_freq_pre.sort_values(ascending=False)
    axs[0].plot(sorted_counts1.index, sorted_counts1.values, marker='o', color='#29B6F6')
    axs[0].set_title('Clonal size of pre-timepoint', fontsize=14)
    axs[0].set_xlabel('Clonal size', fontsize=12)
    axs[0].set_ylabel('Counts', fontsize=12)
    # axs[0].set_xticks(sorted_counts1.index)
    axs[0].tick_params(axis='x', labelsize=10)

    sorted_counts2 = size_freq_pos.sort_values(ascending=False)
    axs[1].plot(sorted_counts2.index, sorted_counts2.values, marker='o', color='#BA68C8')
    axs[1].set_title('Clonal size of post-timepoint', fontsize=14)
    axs[1].set_xlabel('Clonal size', fontsize=12)
    axs[1].set_ylabel('Counts', fontsize=12)
    # axs[1].set_xticks(sorted_counts2.index)
    axs[1].tick_params(axis='x', labelsize=10)

    plt.tight_layout()
    plt.show()

    fig.savefig(savePath, dpi=300)

    return fig


def plotSimilarityCompare(cross_sim, cross_lin_mat, title, savePath):
    fig = plt.figure(figsize=(3.5, 2.5))

    within_clone = np.array(coo_matrix(np.multiply(cross_sim, cross_lin_mat)).data)
    other_value = np.array(coo_matrix(np.multiply(cross_sim, 1 - cross_lin_mat)).data)
    within_clone = within_clone[within_clone <= 0.99]
    other_value = other_value[other_value <= 0.99]
    # print(np.max(within_clone), np.max(other_value))
    colors = ['#98DF8A', '#DBDB8D'] # '#F9DF91', '#A8CFE8'
    plt.hist(other_value, color=colors[0], fill=colors[0], density=True, log=True,
             bins=20, alpha=0.6, label='Across-clone')
    plt.hist(within_clone, color=colors[1], fill=colors[1], density=True, log=True,
             bins=20, alpha=0.6, label='Within-clone')

    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylabel('Density (log-scale)', fontsize=12)
    plt.xlabel('Transcriptomic similarity', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()

    plt.tight_layout()
    plt.savefig(savePath, dpi=300)

    plt.show()
    del(within_clone, other_value)
    gc.collect()

    return fig


def calculateFractions(scd_obj):
    pre_fractions = np.array([sum(scd_obj.data_pre.obs['cluster'] == str(i)) for i in range(scd_obj.n_clus[0])])
    pre_fractions = pre_fractions / sum(pre_fractions)
    pos_fractions = np.array([sum(scd_obj.data_pos.obs['cluster'] == str(i)) for i in range(scd_obj.n_clus[1])])
    pos_fractions = pos_fractions / sum(pos_fractions)
    return pre_fractions, pos_fractions


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def plotFlowSankey(flow_info, pre_colors, pos_colors, pre_fractions=None, pos_fractions=None,
                   label_position='mid', pre_label=None, pos_label=None, label_size=14, title=None,
                   show_cls=True, figwidth=3.5, figheight=4, saveFig=False, saveName='flow.png'):
    pre_names = flow_info['s'].unique()
    pre_names.sort()
    if pre_fractions is None:
        pre_fractions = np.array([flow_info['s_pm'][flow_info['s'] == item].sum() for item in pre_names])

    pos_names = flow_info['t'].unique()
    pos_names.sort()
    if pos_fractions is None:
        pos_fractions = np.array([flow_info['t_pm'][flow_info['t'] == item].sum() for item in pos_names])

    flow_info['s_ix'] = [np.where(pre_names == x)[0][0] for x in flow_info['s']]
    flow_info['t_ix'] = [np.where(pos_names == x)[0][0] for x in flow_info['t']]

    fig, ax = plt.subplots(figsize=(figwidth, figheight), tight_layout=True, facecolor='none')

    pre_x, pos_x = 0, 1
    if label_position == 'mid':
        pre_x, pos_x = 0, 1
    elif label_position == 'twoside':
        pre_x, pos_x = -0.15, 1.15

    for i in range(len(pre_fractions)):
        bottom = pre_fractions[(i + 1):].sum()
        #         bottom = 1 - pre_fractions[0:(i+1)].sum()
        rectangle = ax.bar(x=[0], height=pre_fractions[i], bottom=bottom, color=pre_colors[i],
                           edgecolor='black', fill=True, linewidth=0.7, width=0.16)
        text_y = rectangle[0].get_height() / 2 + bottom
        if show_cls:
            ax.text(x=pre_x, y=text_y, s=str(pre_names[i]), horizontalalignment='center', verticalalignment='center',
                    fontsize=label_size)
    for i in range(len(pos_fractions)):
        bottom = pos_fractions[(i + 1):].sum()
        #         bottom = 1 - pos_fractions[0:(i+1)].sum()
        rectangle = ax.bar(x=[1], height=pos_fractions[i], bottom=bottom, color=pos_colors[i],
                           edgecolor='black', fill=True, linewidth=0.7, width=0.16)
        text_y = rectangle[0].get_height() / 2 + bottom
        if show_cls:
            ax.text(x=pos_x, y=text_y, s=str(pos_names[i]), horizontalalignment='center', verticalalignment='center',
                    fontsize=label_size)

    if pre_label is not None:
        ax.text(x=0, y=-0.05, s=pre_label, horizontalalignment='center', verticalalignment='center', fontsize=14.5)
    if pos_label is not None:
        ax.text(x=1, y=-0.05, s=pos_label, horizontalalignment='center', verticalalignment='center', fontsize=14.5)

    xs = np.linspace(-5, 5, num=100)
    ys = np.array([sigmoid(x) for x in xs])
    xs = xs / 10 + 0.5
    xs *= 0.83
    xs += 0.085
    #     y_start_record = [1 - pre_fractions[0:ii].sum() for ii in range(len(pre_fractions))]
    #     y_end_record = [1 - pos_fractions[0:ii].sum() for ii in range(len(pos_fractions))]
    y_start_record = [pre_fractions[ii:].sum() for ii in range(len(pre_fractions))]
    y_end_record = [pos_fractions[ii:].sum() for ii in range(len(pos_fractions))]
    y_up_start, y_dw_start = 1, 1
    y_up_end, y_dw_end = 1, 1
    axi = 0
    for si in range(len(pre_fractions)):
        cur_flow_info = flow_info.loc[flow_info['s_ix'] == si, :]
        if cur_flow_info.shape[0] > 0:
            for fi in range(cur_flow_info.shape[0]):
                y_up_start = y_start_record[si]
                y_dw_start = y_up_start - cur_flow_info['s_pm'].iloc[fi]
                y_start_record[si] = y_dw_start

                ti = cur_flow_info['t_ix'].iloc[fi]
                y_up_end = y_end_record[ti]
                y_dw_end = y_up_end - cur_flow_info['t_pm'].iloc[fi]
                y_end_record[ti] = y_dw_end

                y_up_start -= 0.0005
                y_dw_start += 0.001
                y_up_end -= 0.0005
                y_dw_end += 0.001

                ys_up = y_up_start + (y_up_end - y_up_start) * ys
                ys_dw = y_dw_start + (y_dw_end - y_dw_start) * ys

                color_s_t = [pre_colors[si], pos_colors[ti]]
                cmap = LinearSegmentedColormap.from_list('mycmap', [color_s_t[0], color_s_t[1]])
                grad_colors = cmap(np.linspace(0, 1, len(xs) - 1))
                grad_colors = [rgb2hex(color) for color in grad_colors]
                for pi in range(len(xs) - 1):
                    ax.fill_between(xs[pi:(pi + 2)], ys_dw[pi:(pi + 2)], ys_up[pi:(pi + 2)], alpha=0.7,
                                    color=grad_colors[pi], edgecolor=None)

            if pre_fractions[(si + 1):].sum() < y_dw_start - 0.001:
                rectangle = ax.bar(x=[0], height=y_dw_start - 0.001 - pre_fractions[(si + 1):].sum() - 0.01,
                                   bottom=pre_fractions[(si + 1):].sum() + 0.005, color='lightgrey',
                                   edgecolor='grey', fill=True, hatch='//', alpha=0.6, linewidth=0.7, width=0.14)

        elif cur_flow_info.shape[0] == 0:
            y_up_start = y_start_record[si]
            y_dw_start = y_up_start - pre_fractions[si]
            y_start_record[si] = y_dw_start

            y_up_end = y_dw_end
            y_dw_end = y_up_end - 0

            y_up_start -= 0.0005
            y_dw_start += 0.001
            y_up_end -= 0.0005
            y_dw_end = y_up_end

            ys_up = y_up_start + (y_up_end - y_up_start) * ys
            ys_dw = y_dw_start + (y_dw_end - y_dw_start) * ys

            ax.fill_between(xs, ys_dw, ys_up, alpha=0.7,
                            color=pre_colors[si])

            color_s_t = [pre_colors[si], 'white']
            cmap = LinearSegmentedColormap.from_list('mycmap', [color_s_t[0], color_s_t[1]])
            grad_colors = cmap(np.linspace(0, 1, len(xs) - 1))
            grad_colors = [rgb2hex(color) for color in grad_colors]
            for pi in range(len(xs) - 1):
                ax.fill_between(xs[pi:(pi + 2)], ys_dw[pi:(pi + 2)], ys_up[pi:(pi + 2)], alpha=0.7,
                                color=grad_colors[pi], edgecolor=None)

    for ti in range(len(pos_fractions)):
        if pos_fractions[(ti + 1):].sum() < y_end_record[ti] - 0.001:
            rectangle = ax.bar(x=[1], height=y_end_record[ti] - 0.001 - pos_fractions[(ti + 1):].sum() - 0.01,
                               bottom=pos_fractions[(ti + 1):].sum() + 0.005, color='lightgrey',
                               edgecolor='grey', fill=True, hatch='//', alpha=0.6, linewidth=0.7, width=0.14)

    for pos in ['right', 'top', 'bottom', 'left']:
        ax.spines[pos].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.set_ylim(-0.01, 1.01)
    ax.patch.set_alpha(0)

    if title is not None:
        ax.set_title(title)

    plt.show()

    if saveFig:
        fig.savefig(saveName, dpi=300, transparent=True, facecolor='none', edgecolor='none', pad_inches=0.0)

    return (fig)


def getColorMap(items, special_case=None, default_palette=("#43D9FE", "#E78AC3", "#FEC643", "#A6D854",
                                                           "#FE6943", "#E5C494", "#33AEB1", "#FFEC1A",
                                                           "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                                                           "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                                                           "#D5BB67", "#6ACC64", "#D65F5F", "#82C6E2",
                                                           "#DC7EC0", "#4878D0", '#B5A2E0', '#F9B475',
                                                           '#50C7CA', '#CF747A', '#63AFF0', '#8792AF', '#E0CB00')):
    # default_palette = sns.color_palette("muted") + sns.color_palette("Set2")
    custom_palette = {}
    items = items.astype('category')
    # print(items)
    for category in items.cat.categories:
        custom_palette[category] = default_palette[len(custom_palette) % len(default_palette)]
    if special_case is not None:
        for case, color in zip(special_case, ['lightgray', '#C0C0C0']):
            custom_palette[case] = color

    return custom_palette


def fateVectorDiversity(adata, cell_2lin_cls, n_neighbors=5, fate_cls_name="Lineage_fate"):
    cfrs_values = []
    n_samples, n_clusters = cell_2lin_cls.shape
    for i in range(n_samples):
        vector = cell_2lin_cls[i, :]
        if np.sum(vector) > 0:
            count_all = np.sum(vector)
            entropy = -sum((count / count_all) * math.log(count / count_all + 1e-9) for count in vector)
            cfrs_values.append(entropy)

    cell_2lin_cls = cell_2lin_cls[adata.obs[fate_cls_name] != 'Missing', :]
    adata = adata[adata.obs[fate_cls_name] != 'Missing']
    cell_2lin_cls = cell_2lin_cls[adata.obs[fate_cls_name] != 'Uncertain', :]
    adata = adata[adata.obs[fate_cls_name] != 'Uncertain']
    data_points = adata.obsm["X_umap"]
    n_samples = len(data_points)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(data_points)
    _, indices = nbrs.kneighbors(data_points)

    afs_values = []
    np.random.seed(123)
    noise = np.random.normal(0, 1e-10, cell_2lin_cls.shape)
    cell_2lin_cls = cell_2lin_cls + noise

    for i in range(n_samples):
        coor_mat = np.corrcoef(cell_2lin_cls[indices[i], :])
        average_fate_similarity = np.mean(coor_mat[0, 1:])
        afs_values.append(average_fate_similarity)

    cfrs, afs = np.mean(cfrs_values), np.mean(afs_values)

    return cfrs, afs


def neighborFateLabelDiversity(adata, n_neighbors=5, fate_cls_name="Lineage_fate"):
    adata = adata[adata.obs[fate_cls_name] != 'Missing']
    adata = adata[adata.obs[fate_cls_name] != 'Uncertain']
    data_points = adata.obsm["X_umap"]
    # Fate cluster
    labels = adata.obs[fate_cls_name]

    n_samples = len(data_points)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(data_points)
    _, indices = nbrs.kneighbors(data_points)

    ncs_values, ecs_values = [], []

    for i in range(n_samples):
        neighbors = indices[i][1:]
        same_label_count = sum(1 for j in neighbors if labels[j] == labels[i])
        ncs_values.append(same_label_count / len(neighbors))

        label_counts = Counter(labels[j] for j in neighbors)
        entropy = -sum((count / len(neighbors)) * math.log(count / len(neighbors) + 1e-9)
                       for count in label_counts.values())
        ecs_values.append(entropy)

    ncs, ecs = np.mean(ncs_values), np.mean(ecs_values)

    return ncs, ecs


def generateDEGs(cur_expr, index, cluster, fate_str, filter_signif = True):
    de_df = pd.DataFrame({'Gene': pd.DataFrame(cur_expr.uns['rank_genes_groups']['names']).iloc[:, index],
                        'Cluster': cluster,
                        'Fate': fate_str, # str(si) + '->' + str(ti) or 'FlowTo' + str(ti)
                        'scores': pd.DataFrame(cur_expr.uns['rank_genes_groups']['scores']).iloc[:, index],
                        'logfoldchanges': pd.DataFrame(cur_expr.uns['rank_genes_groups']['logfoldchanges']).iloc[:, index],
                        'pvals': pd.DataFrame(cur_expr.uns['rank_genes_groups']['pvals']).iloc[:, index],
                        'pvals_adj': pd.DataFrame(cur_expr.uns['rank_genes_groups']['pvals_adj']).iloc[:, index]})
    if filter_signif:
        de_df = de_df[(de_df['pvals_adj'] < 0.05) & (de_df['logfoldchanges'] > 0)]
    else:
        de_df = de_df[de_df['logfoldchanges'] > 0]
    return de_df