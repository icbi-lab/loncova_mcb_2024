#!/usr/bin/env python
# coding: utf-8

"""
Spatial cell-cell interaction analysis using Voronoi diagrams.
Use random permutation to calculate z-scores of interaction
frequencies. Make union over samples

Requirements:
    * Python >= 3.7.0
    * scipy >= 1.4.1
    * seaborn >= 0.10.1
    * pandas >= 1.0.3
    * matplotlib >= 3.2.1
    * numpy >= 1.19.0


Copyright (c) 2024 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""


# Import libraries
#

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mplc
import seaborn as sns
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


sns.set()
import celltypedefs as ctd


# helper function for finding data files
def get_sample_data(dir_name, f_pattern):
    # Get the list of all files in directory tree at given path
    sample_data_files = []
    list_of_files = list()
    for (dir_path, dir_names, file_names) in os.walk(dir_name):
        list_of_files += [os.path.join(dir_path, file) for file in file_names]

    for elem in list_of_files:
        data_file = elem.split("/")[-1]
        if data_file.endswith(f_pattern):
            sample_data_files.append(elem)

    return sample_data_files


def prepare_data_matrix(data_file, remove_phenotype_idx, direction, z_cutoff, one_hot=False):
    sample_pr = data_file.replace("zscores", "pvalues_rightSided")
    sample_pl = data_file.replace("zscores", "pvalues_leftSided")

    sample = pd.read_csv(data_file, sep="\t", index_col=0).to_numpy(dtype="float32")
    sample_pr = pd.read_csv(sample_pr, sep="\t", index_col=0).to_numpy(dtype="float32")
    sample_pl = pd.read_csv(sample_pl, sep="\t", index_col=0).to_numpy(dtype="float32")

    # make mask to filter z-scores for significant enrichment (p < 0.01)
    mask_table = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="float64")

    for i in range(len(ctdef.cellTypeMap.keys())):
        for j in range(len(ctdef.cellTypeMap.keys())):
            if np.abs(sample_pr[i, j]) < 0.01 or np.abs(sample_pl[i, j]) < 0.01:
                mask_table[i, j] = 1

    # copy z-score table and filter with mask
    sample_copy = sample.copy()

    # remove phenotypes we do not want to show in plot

    if len(remove_phenotype_idx) > 0:
        mask_table = np.delete(mask_table, remove_phenotype_idx, axis=0)
        mask_table = np.delete(mask_table, remove_phenotype_idx, axis=1)

        sample_copy = np.delete(sample_copy, remove_phenotype_idx, axis=0)
        sample_copy = np.delete(sample_copy, remove_phenotype_idx, axis=1)

    # and filter with mask
    sample_copy *= mask_table

    sample_copy[np.isnan(sample_copy)] = 0

    if direction == "attract":
        sample_copy[sample_copy < z_cutoff] = 0
    elif direction == "avoid":
        sample_copy[sample_copy > z_cutoff * -1] = 0
    elif direction == "both":
        sample_copy[abs(sample_copy) < z_cutoff] = 0
    else:
        print("Error no such direction: " + direction)
        sys.exit(1)

    if one_hot:
        sample_copy[sample_copy >= z_cutoff] = 1
        sample_copy[sample_copy <= z_cutoff * -1] = -1
        sample_copy[np.isposinf(sample_copy)] = 1
        sample_copy[np.isneginf(sample_copy)] = -1
    else:
        # deal with inf
        sample_copy[np.isposinf(sample_copy)] = np.nan
        sample_copy[np.isnan(sample_copy)] = np.nanmax(sample_copy) + 0.1
        sample_copy[np.isneginf(sample_copy)] = np.nan
        sample_copy[np.isnan(sample_copy)] = np.nanmin(sample_copy) - 0.1

    return sample_copy


# main
if __name__ == "__main__":

    # Argument parser
    parser = argparse.ArgumentParser(description="Run spatial cell-cell interaction one-hot encoded z-score clustering")
    parser.add_argument("--data_dir", required=True, default="", type=str, help="Data directory")
    parser.add_argument("--result_dir", required=True, default="", type=str, help="Results directory")

    parser.add_argument(
        "--celltype_def", required=True, default="", type=str, help="Celltype name/color/number definitions in tsv file"
    )

    parser.add_argument(
        "--granularity",
        required=False,
        default="full",
        type=str,
        choices=["full", "lumped"],
        help="Set celltype granularity",
    )

    parser.add_argument(
        "--onehot",
        dest="one_hot",
        action="store_true",
        help="Use one hot encoding",
    )
    parser.set_defaults(one_hot=False)

    parser.add_argument(
        "--s_count",
        dest="s_count",
        action="store_true",
        help="Count samples with significant celltype pair attraction/avoidance",
    )
    parser.set_defaults(s_count=False)

    parser.add_argument(
        "--direction",
        choices=["attract", "avoid", "both"],
        type=str,
        default="both",
        help="Interaction direction",
    )

    parser.add_argument(
        "--max_z",
        required=False,
        default=130,
        type=int,
        help="Maximum abs(z-score) value to anchor the colormap (default 130)",
    )

    parser.add_argument(
        "--z_cutoff",
        required=False,
        default=2,
        type=int,
        help="Minimum abs(z-score) value to consider as attraction or avoidance (default 2)",
    )

    parser.add_argument(
        "--s_count_min",
        required=False,
        default=3,
        type=int,
        help="Minimum number of samples that attraction/avoidance should be labelled",
    )

    parser.add_argument(
        "--discrete_cm",
        dest="discrete_cm",
        action="store_true",
        help="Use discrete color map with 15 categories in heatmap plots",
    )
    parser.set_defaults(discrete_cm=False)

    parser.add_argument(
        "--num_categories",
        required=False,
        default=10,
        type=int,
        help="Number of categories (one direction) for discrete colormap (default 10)",
    )

    parser.add_argument(
        "--legend",
        dest="legend",
        action="store_true",
        help="Add phenotype legend (default false)",
    )
    parser.set_defaults(legend=False)

    parser.add_argument(
        "--legend_loc",
        required=False,
        choices=["bottom", "right"],
        type=str,
        default="bottom",
        help="Phenotype legend location (default bottom)",
    )

    parser.add_argument(
        "--mark_celltype_classes",
        dest="mark_celltype_classes",
        action="store_true",
        help="Mark cell type classes with class color (default false)",
    )
    parser.set_defaults(mark_celltype_classes=False)

    parser.add_argument(
        "--skip_phenotype",
        required=False,
        default="",
        type=str,
        nargs="*",
        help="Do not plot listed phenotypes",
    )

    parser.add_argument(
        "--skip_excluded", required=False, action="store_true", help="skip plotting phenotypes with cname excluded"
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    # Input output settings
    result_dir = args.result_dir
    data_dir = args.data_dir
    celltype_def_file = args.celltype_def
    granularity = args.granularity

    # initialize the celltype defs
    ctdef = ctd.ctdef(granularity, celltype_def_file)

    # set direction of interaction
    direction = args.direction

    # use one hot encoding
    one_hot = args.one_hot

    # for each cell type pair count samples with significant attraction/avoidance
    s_count = args.s_count
    # we need ones for each significant attraction/avoidance if counting samples
    one_hot = s_count

    # Min/Max z-score values to consider and anchor the colormap
    vmin = args.max_z * -1
    vmax = args.max_z

    # z-score cutoff
    z_cutoff = args.z_cutoff

    # min number of samples needed to label attraction or avoidance
    s_count_min = args.s_count_min

    # use discrete cm
    discrete_cm = args.discrete_cm
    num_categories = args.num_categories

    if direction == "both":
        title_direction = "attraction/avoidance"
    elif direction == "attract":
        title_direction = "attraction"
    elif direction == "avoid":
        title_direction = "avoidance"

    result_name = "summed_sample_zscores_" + title_direction.replace("/", "_")
    if s_count:
        result_name = "sample_count_" + title_direction.replace("/", "_")

    # legend
    add_legend = args.legend
    legend_loc = args.legend_loc

    # cell type classes
    mark_celltype_classes = args.mark_celltype_classes

    # list of phenotypes to skip in plots
    skip_phenotype = args.skip_phenotype
    skip_excluded = args.skip_excluded
    skipped_phenotypes = False

    #
    # Start processing data files
    #

    z_score_data_processed = []

    z_score_datafiles = get_sample_data(data_dir, "_voronoi_neighbors_zscores.tsv")
    for data_file in z_score_datafiles:
        print(data_file)

        # remove phenotypes we do not want to show in plot
        if len(skip_phenotype) > 0 or skip_excluded:
            remove_phenotype_idx = []

            if skip_phenotype == "":
                skip_phenotype = []

            for ct in ctdef.cellTypeNameMap.keys():
                if ctdef.cellTypeExclude[ct]:
                    skip_phenotype.append(ct)

            skip_phenotype = list(set(skip_phenotype))

            for phenotype in skip_phenotype:
                if phenotype in ctdef.cellTypeMap:
                    remove_phenotype_idx.append(ctdef.cellTypeMap[phenotype])

                    # only do this for the first image
                    if skipped_phenotypes is not True:
                        color_map_idx = ctdef.cellTypeNameColorMapIdx[ctdef.cellTypeNameMap[phenotype]]
                        if color_map_idx in ctdef.cellTypeIdxColorMap:
                            del ctdef.cellTypeIdxColorMap[color_map_idx]
                            del ctdef.cellTypeColorMap[phenotype]
                            ctdef.cellTypeMapOrder_idx.remove(color_map_idx)

            # done skipping phenotypes
            skipped_phenotypes = True

        tmp_data = prepare_data_matrix(data_file, remove_phenotype_idx, direction, z_cutoff, one_hot)
        z_score_data_processed.append(tmp_data)

    sample_sum = np.nansum(z_score_data_processed, axis=0)

    # Do Stouffer
    if not one_hot:
        sample_sum = sample_sum / len(z_score_datafiles) ** 0.5

    # set vmin/vmax to number of samples if counting samples
    if s_count:
        vmax = len(z_score_datafiles)
        vmin = vmax * -1
        if discrete_cm:
            num_categories = len(z_score_datafiles)
        if direction == "avoid":
            sample_sum = sample_sum * -1

    # assign cell type number to rows and columns as used in manuscript
    sample_sum_df = pd.DataFrame(
        sample_sum, columns=ctdef.cellTypeIdxColorMap.keys(), index=ctdef.cellTypeIdxColorMap.keys()
    )

    plot_title = "Cell type interactions (" + title_direction + "): all patients summed"
    if s_count:
        plot_title = "Significant cell type interactions (" + title_direction + "): patient counts"

    # add legend
    if add_legend:
        # set cell type names and colors for legend
        patch_list = []
        for key, val in ctdef.cellTypeColorMap.items():
            data_key = mpatches.Patch(
                color=val,
                label=str(
                    ctdef.cellTypeNameColorMapIdx[ctdef.cellTypeNameMap[key]] + ": " + ctdef.cellTypeNameMap[key]
                ),
            )
            patch_list.append(data_key)

    #
    # plot the custom ordered heatmap using seaborn
    #

    # color map
    cbar_kws = {"shrink": 0.7}
    cbar_pos = (1.0, 0.15, 0.03, 0.5)

    if direction == "both":
        cmap = plt.get_cmap("seismic")
    elif direction == "attract":
        vmin = z_cutoff
        cmap = plt.get_cmap("Reds")
        cmap.set_under("white")
    elif direction == "avoid" and not s_count:
        vmin = vmax * -1
        vmax = z_cutoff * -1
        cmap = plt.get_cmap("Blues_r")
        cmap.set_over("white")
    elif direction == "avoid" and s_count:
        vmin = 0
        cmap = plt.get_cmap("Blues")
        cmap.set_under("white")
    else:
        print("Ouch")
        sys.exit(1)

    # discretized with 15 catergories
    if discrete_cm:
        if direction == "both":
            num_categories = num_categories * 2 + 1
            c = sns.color_palette("seismic", num_categories)
        elif direction == "attract" and not s_count:
            num_categories = num_categories + 1
            c = sns.color_palette("Reds", num_categories)
        elif direction == "attract" and s_count:
            vmin = 1
            num_categories = num_categories
            c = sns.color_palette("Reds", num_categories)
            cmap.set_under("white")
        elif direction == "avoid" and not s_count:
            num_categories = num_categories + 1
            c = sns.color_palette("Blues_r", num_categories)
        elif direction == "avoid" and s_count:
            vmin = 1
            num_categories = num_categories
            c = sns.color_palette("Blues", num_categories)
            cmap.set_under("white")
        else:
            print("Ouch")
            sys.exit(1)

        cmap = mplc.LinearSegmentedColormap.from_list("My_discrete", c, num_categories)
        if direction != "both":
            cmap.set_under("white")

    if s_count:
        annotation = np.where(abs(sample_sum) >= s_count_min, abs(sample_sum).astype(int).astype(str), "")
    else:
        annotation = False
        a_mask = False

    if one_hot and not s_count:
        if direction == "attract":
            cmap = [(1, 1, 1), (1, 0, 0)]
            vmax = 1
            vmin = 0
            cbar_kws = {"shrink": 0.2, "ticks": [0, 1]}
            cbar_pos = (1.0, 0.4, 0.03, 0.05)
        elif direction == "avoid":
            cmap = [(0, 0, 1), (1, 1, 1)]
            vmax = 0
            vmin = -1
            cbar_kws = {"shrink": 0.2, "ticks": [-1, 0]}
            cbar_pos = (1.0, 0.4, 0.03, 0.05)
        elif direction == "both":
            cmap = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]
            vmax = 1
            vmin = -1
            cbar_kws = {"shrink": 0.2, "ticks": [-1, 0, 1]}
            cbar_pos = (1.0, 0.4, 0.03, 0.05)

    # with sns.axes_style("white"):
    fig = sns.clustermap(
        sample_sum_df,
        method="weighted",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        cbar_pos=cbar_pos,
        cbar_kws=cbar_kws,
        annot=annotation,
        fmt="",
        yticklabels=True,
        xticklabels=True,
        row_colors=ctdef.cellTypeIdxColorMap_s,
        col_colors=ctdef.cellTypeIdxColorMap_s,
        linewidths=0.5,
        # linecolor=[0.85, 0.85, 0.85],
        figsize=(12, 12),
    )

    if one_hot and not s_count:
        if direction == "both":
            fig.ax_cbar.set_yticklabels(["avoidance", "", "attraction"])
        elif direction == "attract":
            fig.ax_cbar.set_yticklabels(["", "attraction"])
        elif direction == "avoid":
            fig.ax_cbar.set_yticklabels(["avoidance", ""])
    elif s_count:
        fig.ax_cbar.set_label("# samples")
        labels = np.array([["A", "B"], ["C", "D"], ["E", "F"]])
    else:
        fig.ax_cbar.set_label("z-score")

    fig.ax_heatmap.tick_params(axis="both", which="both", length=0, pad=28, labelsize="x-small")
    fig.ax_heatmap.yaxis.set_ticks_position("left")
    fig.ax_heatmap.xaxis.set_ticks_position("top")

    fig.ax_col_colors.set_yticks([])
    fig.ax_row_colors.set_xticks([])

    fig.ax_row_dendrogram.set_visible(False)
    fig.ax_col_dendrogram.set_visible(False)

    fig.ax_heatmap.set_title(plot_title, pad=20)

    # add legend
    if add_legend:
        if legend_loc == "bottom":
            l_loc = "upper center"
            ncol = 5
            bbox_to_anchor = (0.5, -0.05)
        if legend_loc == "right":
            l_loc = "center right"
            ncol = 1
            bbx = 1.5 if mark_celltype_classes else 1.45
            bbox_to_anchor = (bbx, 0.5)

        fig.ax_heatmap.legend(
            handles=patch_list,
            loc=l_loc,
            bbox_to_anchor=bbox_to_anchor,
            fontsize="xx-small",
            ncol=ncol,
            labelspacing=0.25,
            handletextpad=0.5,
            facecolor="white",
            frameon=False,
        )

    plt.tight_layout()

    if one_hot and not s_count:
        result_name += "_oneHot"

    fig.savefig(
        result_dir + "/" + result_name + "_clustered.png",
        dpi=300,
    )

    #
    # plot the custom ordered heatmap using seaborn
    #

    # order z-score table by custom order
    custom_sorted_sample_sum_df = sample_sum_df[ctdef.cellTypeMapOrder_idx]
    custom_sorted_sample_sum_df = custom_sorted_sample_sum_df.reindex(ctdef.cellTypeMapOrder_idx)

    # compose row/col colors and labels
    if mark_celltype_classes:
        row_colors = pd.concat(
            [ctdef.cellTypeDummyIdxColorMap_s, ctdef.cellTypeIdxColorMap_s, ctdef.cellTypeClassIdxColorMap_s], axis=1
        )
    else:
        row_colors = ctdef.cellTypeIdxColorMap_s
    col_colors = row_colors

    fig = sns.clustermap(
        custom_sorted_sample_sum_df,
        row_cluster=False,
        col_cluster=False,
        cmap=cmap,
        vmin=vmin,  # -(np.max(pair_z_score_table_copy))*0.75,
        vmax=vmax,  # np.max(pair_z_score_table_copy)*0.75,
        cbar_pos=cbar_pos,
        cbar_kws=cbar_kws,
        annot=annotation,
        fmt="",
        yticklabels=True,
        xticklabels=True,
        row_colors=row_colors,
        col_colors=col_colors,
        linewidths=0.5,
        figsize=(12, 12),
    )

    if one_hot and not s_count:
        if direction == "both":
            fig.ax_cbar.set_yticklabels(["avoidance", "", "attraction"])
        elif direction == "attract":
            fig.ax_cbar.set_yticklabels(["", "attraction"])
        elif direction == "avoid":
            fig.ax_cbar.set_yticklabels(["avoidance", ""])
    elif s_count:
        fig.ax_cbar.set_label("# samples")
        labels = np.array([["A", "B"], ["C", "D"], ["E", "F"]])
    else:
        fig.ax_cbar.set_label("z-score")

    if mark_celltype_classes:
        fig.ax_heatmap.tick_params(axis="both", which="both", length=0, pad=56, labelsize="x-small")
    else:
        fig.ax_heatmap.tick_params(axis="both", which="both", length=0, pad=28, labelsize="x-small")

    fig.ax_heatmap.yaxis.set_ticks_position("left")
    fig.ax_heatmap.xaxis.set_ticks_position("top")

    fig.ax_col_colors.set_yticks([])
    fig.ax_row_colors.set_xticks([])

    fig.ax_row_dendrogram.set_visible(False)
    fig.ax_col_dendrogram.set_visible(False)

    if mark_celltype_classes:
        fig.ax_heatmap.set_title(plot_title, pad=40)
    else:
        fig.ax_heatmap.set_title(plot_title, pad=20)

    # add legend
    if add_legend:
        if legend_loc == "bottom":
            l_loc = "upper center"
            ncol = 5
            bbox_to_anchor = (0.5, -0.05)
        if legend_loc == "right":
            l_loc = "center right"
            ncol = 1
            bbx = 1.5 if mark_celltype_classes else 1.45
            bbox_to_anchor = (bbx, 0.5)

        fig.ax_heatmap.legend(
            handles=patch_list,
            loc=l_loc,
            bbox_to_anchor=bbox_to_anchor,
            fontsize="xx-small",
            ncol=ncol,
            labelspacing=0.25,
            handletextpad=0.5,
            facecolor="white",
            frameon=False,
        )

    plt.tight_layout()

    if one_hot and not s_count:
        result_name += "_oneHot"

    fig.savefig(
        result_dir + "/" + result_name + "_ordered.png",
        dpi=300,
    )

    print(result_name)
