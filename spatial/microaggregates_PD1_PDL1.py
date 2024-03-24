#!/usr/bin/env python
# coding: utf-8

"""
Analyze microaggregates via Voronoi diagrams of phenotyped IMC images.

Requirements:
    * Python >= 3.7.0
    * scipy >= 1.4.1
    * pandas == 1.3.5
    * matplotlib >= 3.2.1
    * numpy >= 1.19.0
    * networkx


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

# Import python packages


import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.patches as mpatches
from scipy.spatial import Voronoi
from scipy import stats
import networkx as nx
from multiprocessing import Pool

from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

import celltypedefs as ctd

# sample assignment
def get_sample_map(sample_sheet):
    sample_map = {}
    sample_sheet = pd.read_csv(sample_sheet, sep="\t")

    sample_map = sample_sheet.set_index('sample_ID')['sample_name'].to_dict()

    return sample_map


# helper function for finding data files


def get_sample_data(dir_name, sample_id):
    # Get the list of all files in directory tree at given path
    sample_data_files = []
    list_of_files = list()
    for (dir_path, dir_names, file_names) in os.walk(dir_name):
        list_of_files += [os.path.join(dir_path, file) for file in file_names]

    for elem in list_of_files:
        data_file = elem.split("/")[-1]
        if data_file.startswith(sample_id):
            sample_data_files.append(elem)

    return sample_data_files


# Function to identify microaggregates from neighboring cells


def get_microaggregates(contacts_oi, phenotypes):
    # Make a graph from the touching cell pairs
    G = nx.Graph()
    G.add_edges_from(contacts_oi)

    # Get a list of connected components (list of lists). Each list is a connected cell cluster aka microaggregate
    connected_components = list(nx.connected_components(G))

    # Filter connected cell clusters with less than `min_members` cells to get microaggregates
    microaggregates = [sub for sub in connected_components if len(sub) >= min_members]

    # Remove microaggregates that consist only of a single cell category
    microaggregates_final = []

    for sub in microaggregates:
        tc_count = 0
        ic_count = 0
        for cell in sub:
            tc_count += 1 if phenotypes["tumor_pheno"][cell] == 1 else 0
            ic_count += 1 if phenotypes["immune_pheno"][cell] == 1 else 0
        if tc_count != 0 and ic_count != 0:
            microaggregates_final.append(list(sub))

    return microaggregates_final


def identify_microaggregates(vor, coordinates, rnd_phenotypes=None):
    contacts_oi = []

    if rnd_phenotypes is not None:
        phenotypes = rnd_phenotypes
    else:
        phenotypes = coordinates[["phenotype", "category", "PD1_marker", "PDL1_marker", "immune_pheno", "tumor_pheno"]]

    for interaction in vor.ridge_points:
        cell_a = interaction[0]
        cell_b = interaction[1]

        if cell_a <= max_cell_idx and cell_b <= max_cell_idx:

            cell_a_pos = np.asarray(
                [
                    coordinates["Location_Center_X"][cell_a],
                    coordinates["Location_Center_Y"][cell_a],
                ]
            )
            cell_b_pos = np.asarray(
                [
                    coordinates["Location_Center_X"][cell_b],
                    coordinates["Location_Center_Y"][cell_b],
                ]
            )

            cell_dist = np.linalg.norm(cell_a_pos - cell_b_pos)

            if cell_dist < max_dist:
                if strict:
                    if (
                        phenotypes["immune_pheno"][cell_a] == 1
                        and phenotypes["PD1_marker"][cell_a] == 1
                        and phenotypes["PDL1_marker"][cell_a] == 0
                        and phenotypes["tumor_pheno"][cell_b] == 1
                        and phenotypes["PDL1_marker"][cell_b] == 1
                        and phenotypes["PD1_marker"][cell_b] == 0
                    ) or (
                        phenotypes["immune_pheno"][cell_b] == 1
                        and phenotypes["PD1_marker"][cell_b] == 1
                        and phenotypes["PDL1_marker"][cell_b] == 0
                        and phenotypes["tumor_pheno"][cell_a] == 1
                        and phenotypes["PDL1_marker"][cell_a] == 1
                        and phenotypes["PD1_marker"][cell_a] == 0
                    ):
                        contacts_oi.append(interaction)
                else:
                    if (
                        phenotypes["immune_pheno"][cell_a] == 1
                        and phenotypes["PD1_marker"][cell_a] == 1
                        and phenotypes["tumor_pheno"][cell_b] == 1
                        and phenotypes["PDL1_marker"][cell_b] == 1
                    ) or (
                        phenotypes["immune_pheno"][cell_b] == 1
                        and phenotypes["PD1_marker"][cell_b] == 1
                        and phenotypes["tumor_pheno"][cell_a] == 1
                        and phenotypes["PDL1_marker"][cell_a] == 1
                    ):
                        contacts_oi.append(interaction)

    contacts_oi = np.asarray(contacts_oi)
    microaggregates_final = get_microaggregates(contacts_oi, phenotypes)

    return microaggregates_final, len(contacts_oi)


# Function to generate random PDL1+ tumor / PD1+ immune microaggregates
# We keep cell numbers, cell types and spatial organization as in the original image
# but randomly re-assign the cell types (pheontype)


def register_rnd_microaggregates(p):
    # random sample phenotypes
    rng = np.random.RandomState()
    rnd_phenotypes = (
        coordinates[["phenotype", "category", "PD1_marker", "PDL1_marker", "immune_pheno", "tumor_pheno"]]
        .sample(frac=1, random_state=rng)
        .reset_index(drop=True)
    )

    rnd_microaggregates_final, rnd_contacts_oi_counts = identify_microaggregates(vor, coordinates, rnd_phenotypes)

    return rnd_microaggregates_final, rnd_contacts_oi_counts


# ###   function to calculate the z-score
#
# for p-value calculation we use:
#
# Phipson, B. & Smyth, G. K. Permutation P-values should never be zero:
# calculating exact P-values when permutations are randomly drawn. Stat.
# Appl. Genet. Mol. Biol. 9, Article39 (2010).
#
def calculate_zscore_stat(sum_close_pairs, rnd_close_pairs_counts, num_permutations):
    z = (sum_close_pairs - np.nanmean(rnd_close_pairs_counts)) / np.nanstd(rnd_close_pairs_counts)
    p = stats.norm.sf(abs(z))

    rnd_greater = 0
    if not np.isnan(sum_close_pairs):
        rnd_greater = np.sum([i > sum_close_pairs for i in rnd_close_pairs_counts])
        rnd_smaller = np.sum([i < sum_close_pairs for i in rnd_close_pairs_counts])
        pr = (rnd_greater + 1) / (num_permutations + 1)
        pl = (rnd_smaller + 1) / (num_permutations + 1)

    return z, p, pr, pl



if __name__ == "__main__":

    # Argument parser
    parser = argparse.ArgumentParser(description="Analyze PD1/PDL1 microaggregates with Voronoi diagrams")
    parser.add_argument("--result_dir", required=True, default="", type=str, help="Results directory")
    parser.add_argument(
        "--data_dir",
        required=True,
        default="",
        type=str,
        help='Data directory with phenotype coordinate files: i.e. path to "Single_cell_data" directory',
    )
    parser.add_argument(
        "--run_samples",
        required=False,
        default=-1,
        type=int,
        help="Run only this many samples (default: all)",
    )
    parser.add_argument(
        "--run_sample",
        required=False,
        default="",
        type=str,
        help="Run only this samples (default: all)",
    )
    parser.add_argument(
        "--sample_sheet",
        required=True,
        default="",
        type=str,
        help='TSV file with two columns <sample_ID> and <sample_name>: maps sample IDs to sample names shown in results',
    )
    parser.add_argument(
        "--celltype_def", required=True, default="", type=str, help="Celltype name/color/number definitions in tsv file"
    )
    parser.add_argument(
        "--min_agg_size",
        required=False,
        default=2,
        type=int,
        help="Minimum cells in an microaggregate (default: 2)",
    )
    parser.add_argument(
        "--strict",
        required=False,
        default=False,
        action="store_true",
        help="Consider only strict pairs of both immune PD1+PDL1- / tumor PD1-PDL1+ in an aggregate (default: True)",
    )
    parser.add_argument(
        "--max_dist",
        required=False,
        default=25,
        type=int,
        help="Maximum center to center distance between two cells in an microaggregate (default: 25)",
    )
    parser.add_argument(
        "--plot_size_x",
        required=False,
        default=10,
        type=float,
        help="Size of plot in inches horizontal",
    )
    parser.add_argument(
        "--plot_size_y",
        required=False,
        default=10,
        type=float,
        help="Size of plot in inches vertical",
    )
    parser.add_argument(
        "--image_size_x",
        required=False,
        default=1000,
        type=int,
        help="Size of original ROI image in pixels horizontal",
    )
    parser.add_argument(
        "--image_size_y",
        required=False,
        default=1000,
        type=int,
        help="Size of original ROI image in pixels vertical",
    )
    parser.add_argument(
        "--make_cropped_plot",
        required=False,
        default=False,
        action="store_true",
        help="Make a cropped region of the plot, it set default origin is 0,0 and copped size is 200x200, may be adjusted via --image_crop_origin x,y and image_crop_size x,y",
    )
    parser.add_argument(
        "--image_crop_origin",
        required=False,
        default=[0, 0],
        nargs=2,
        type=int,
        help="x y coordinates of origin for cropped image",
    )
    parser.add_argument(
        "--image_crop_size",
        required=False,
        default=[200, 200],
        nargs=2,
        type=int,
        help="x y size in pixels of cropped image",
    )
    parser.add_argument(
        "--ridge_line_width",
        required=False,
        default=0.25,
        type=float,
        help="line width of Voronoi ridges",
    )
    parser.add_argument(
        "--nucleus_marker_size",
        required=False,
        default=0.5,
        type=float,
        help="size of Voronoi point marker",
    )
    parser.add_argument(
        "--num_cores",
        required=False,
        default=1,
        type=int,
        help="number of cores to use for parallel computing",
    )
    parser.add_argument(
        "--num_permutations",
        required=False,
        default=1000,
        type=int,
        help="Number of random permutations to run",
    )
    parser.add_argument(
        "--no_annotation",
        required=False,
        default=False,
        action="store_true",
        help="Do not add legend or title",
    )
    parser.add_argument(
        "--fig_format",
        required=False,
        choices=["svg", "pdf", "png"],
        default="png",
        help="Do not add legend or title",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    # Minimum number of members of an aggregates
    min_members = args.min_agg_size

    # Use strict aggregation mode: only cell paris consisting of PD1+PDL1- iummune cells and PD1-PDL1+ tumor cells are considered
    strict = args.strict

    # Maximum distance between cell pairs
    max_dist = args.max_dist

    # Input output settings
    result_dir = args.result_dir
    data_dir = args.data_dir
    sample_sheet = args.sample_sheet
    celltype_def_file = args.celltype_def

    # plot in inches
    x_max = args.plot_size_x
    y_max = args.plot_size_y

    # image size in px
    x_max_img = args.image_size_x
    y_max_img = args.image_size_y

    # set voronoi ridge line width
    ridge_line_width = args.ridge_line_width

    # set voronoi point size
    nucleus_marker_size = args.nucleus_marker_size

    # number of cores to use in parallel functions
    num_cores = args.num_cores

    # number of random permutations
    num_permutations = args.num_permutations

    # lets see if we need to make a cropped plot
    if args.make_cropped_plot:
        image_crop_origin = args.image_crop_origin
        image_crop_size = args.image_crop_size

        if (
            image_crop_origin[0] + image_crop_size[0] >= x_max_img
            or image_crop_origin[1] + image_crop_size[1] >= y_max_img
        ):
            print("Cropped section must fit into ROI image")
            sys.exit(1)

        if image_crop_size[0] != image_crop_size[1]:
            if image_crop_size[0] > image_crop_size[1]:
                scale_plot_y = image_crop_size[1] / image_crop_size[0]
                y_max = y_max * scale_plot_y
            else:
                scale_plot_x = image_crop_size[0] / image_crop_size[1]
                x_max = x_max * scale_plot_x

        ridge_line_width = 1000 / max(image_crop_size) * ridge_line_width
        nucleus_marker_size = 1000 / max(image_crop_size) * nucleus_marker_size

    # stop after running this many samples (-1 runs all)
    nr_of_samples_to_run = args.run_samples

    # only run this sample
    run_sample = args.run_sample

    # annotation
    no_annotation = args.no_annotation

    # figure store format
    fig_format = args.fig_format

    # image/sample counter
    i_counter = 0

    # get the sample_ID to sample_name mapping
    # sample ids are the ids/names that the filenames with the
    # cell coordinates for a specific sample are starting with
    # one sample can have mutiple coordinates files (ROIs)
    # sample names are used for all plots and tables.
    sample_map = get_sample_map(sample_sheet)

    # List for final aggregate sizes and counts
    aggregate_sizes = []
    aggregate_counts = []
    pair_oi_counts = []

    # List for PD1, PDL1 cell counts
    cell_counts = []

    # initialize the celltype defs
    ctdef = ctd.ctdef("full", celltype_def_file)

    pb_samples = tqdm(sample_map.keys())
    for s in pb_samples:
        pb_samples.set_description("Processing %s sample" % len(sample_map))

        # only run specified sample
        if run_sample != "" and sample_map[s] != run_sample:
            continue

        # stop after nr_of_samples_to_run (if requested)
        if i_counter == nr_of_samples_to_run and nr_of_samples_to_run > 0:
            print("Done " + str(i_counter) + " images.")
            break

        # reset count tables
        sample_microaggregate_table = []
        sample_rnd_microaggregate_count_table = []
        sample_contacts_oi_count_table = []
        sample_rnd_contacts_oi_count_table = []

        sample_id = sample_map[s]
        # get images for sample
        sample_data_files = get_sample_data(data_dir, s)

        # iterate over sample ROI images
        pb_sample_data_files = tqdm(sample_data_files)
        for sample_data_file in pb_sample_data_files:
            fname = os.path.basename(sample_data_file)
            pb_sample_data_files.set_description("Processing %s - %s" % (sample_map[s], fname))

            tile_name = fname.replace("__PD1_PDL1.csv", "")

            # read ROI coordinate file
            coordinates = pd.read_csv(sample_data_file, na_values="NA", low_memory=False)

            # get maximum index of cells
            max_cell_idx = len(coordinates) - 1

            # get cell positions in ROI to calculate Voronoi diagram
            positions = np.transpose(
                np.array(
                    [
                        coordinates["Location_Center_X"].values,
                        coordinates["Location_Center_Y"].values,
                    ]
                )
            )

            # extend the positions beyond the bounding box of the image
            positions = np.append(
                positions,
                [[9999, 9999], [-9999, 9999], [9999, -9999], [-9999, -9999]],
                axis=0,
            )

            # calculate Voronoi
            vor = Voronoi(positions)

            # Get all contacts of interest from Voronoi ridge points. Cells of a touching cell pair have to be tumor_pheno PDL1+
            # immune_pheno PD1+, if `strict == True` we only accept pairs in which are tumor PDL1+PD1- and immune PD1+PDL1-
            # max_dist is maximum distance between center of mass in a pair (diameter is typically 10-25 um)
            microaggregates_final, contacts_oi_counts = identify_microaggregates(vor, coordinates, rnd_phenotypes=None)

            # append microaggregates of ROI to sample microaggregate_table
            sample_microaggregate_table.append(microaggregates_final)

            # append contacts of interest of ROI to sample_contacts_oi_count_table
            sample_contacts_oi_count_table.append(contacts_oi_counts)

            # set up process pool for parallel processing and run Monte Carlo
            pool = Pool(num_cores)
            rnd_microaggregates, rnd_contacts_oi_counts = zip(
                *pool.map(register_rnd_microaggregates, [p for p in range(num_permutations)])
            )
            pool.close()

            # store monte carlo results for current image ROI
            # rnd_microaggregates = np.asarray(rnd_microaggregates)
            rnd_microaggregates = np.sum(np.array(rnd_microaggregates, dtype=object), axis=0)
            rnd_microaggregate_counts = [len(i) for i in rnd_microaggregates]

            # append monte carlo results of current image ROI to sample monte carlo results
            sample_rnd_microaggregate_count_table.append(rnd_microaggregate_counts)

            # append rnd contacts of interest of ROI to sample_contacts_oi_count_table
            sample_rnd_contacts_oi_count_table.append(rnd_contacts_oi_counts)

            ### Make Voronoi plot with cells in microaggregates marked by cell type color

            # get list of cells that are member of a filtered microaggregate (flattens the list of lists)
            microaggreagate_cells = sum(microaggregates_final, [])

            # make Voronoi diagram figure
            p = plt.figure(figsize=(x_max, y_max))

            # define cell nucleus center of mass marker
            marker_style = dict(
                color="tab:blue",
                marker="o",
                markersize=nucleus_marker_size,
                fillstyle="full",
            )

            # reset
            cell_color = {}
            cell_nr = 0
            cell_types_present = []
            patch_list = []

            # loop over point_regions and draw corresponding polygons from vertices
            # and fill with cell type specific color
            for p_reg in vor.point_region:
                if not -1 in vor.regions[p_reg] and cell_nr <= max_cell_idx:
                    polygon = [vor.vertices[i] for i in vor.regions[p_reg]]

                    # first fill with unknown cell type color
                    plt.fill(
                        # *zip(*polygon), ctdef.cellTypeColorMap["Unknown"], edgecolor="w", linewidth=ridge_line_width
                        *zip(*polygon),
                        "white",
                        edgecolor="r",
                        linewidth=ridge_line_width,
                    )

                    # if cell is in microaggregate color it by the cell type specific color
                    if cell_nr in microaggreagate_cells:
                        cell_type = coordinates["phenotype"][cell_nr]

                        c_col = ctdef.cellTypeColorMap[cell_type]
                        plt.fill(
                            *zip(*polygon),
                            c_col,
                            edgecolor="w",
                            linewidth=ridge_line_width,
                        )
                        cell_types_present.append(cell_type)

                cell_nr += 1

            # plot nuclei center of mass
            plt.plot(positions[:-1, 0], positions[:-1, 1], "o", **marker_style)

            # generate celltype legend
            cell_types_present = set(cell_types_present)

            for cell_type in cell_types_present:
                data_key = mpatches.Patch(
                    color=ctdef.cellTypeColorMap[cell_type],
                    label=ctdef.cellTypeNameMap[cell_type],
                )
                patch_list.append(data_key)

            # add title
            if not no_annotation:
                plt.title(sample_id + ": Tumor PDL1+ / Immune PD1+ microaggregates")  # title with fontsize 20

            plt.xlim([0, x_max_img]), plt.ylim([0, y_max_img])

            if args.make_cropped_plot:
                # plot only cropped section
                plt.xlim([image_crop_origin[0], image_crop_origin[0] + image_crop_size[0]])
                plt.ylim([image_crop_origin[1], image_crop_origin[1] + image_crop_size[1]])
                plt.xticks(np.arange(image_crop_origin[0], image_crop_origin[0] + image_crop_size[0]+1, step=100), np.arange(0,image_crop_size[0]+1, step=100))
                plt.yticks(np.arange(image_crop_origin[1], image_crop_origin[1] + image_crop_size[1]+1, step=100), np.arange(0,image_crop_size[1]+1, step=100))

            # add legend
            if not no_annotation:
                plt.legend(
                    handles=patch_list,
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.05),
                    fontsize="xx-small",
                    ncol=5,
                    labelspacing=0.25,
                    handletextpad=0.5,
                )
            else:
                plt.xticks([image_crop_origin[0], image_crop_origin[0] + image_crop_size[0]], [0, image_crop_size[0]])
                plt.yticks([image_crop_origin[1], image_crop_origin[1] + image_crop_size[1]], [0, image_crop_size[1]])

                plt.tick_params(top=False, bottom=False, left=False, right=False,
                                labelleft=True, labelbottom=True)


            # plt.show()

            plt.tight_layout()

            aggregate_mode = "_relaxed"
            if strict:
                aggregate_mode = "_strict"

            if args.make_cropped_plot:
                plt.savefig(
                    result_dir
                    + "/"
                    + sample_id
                    + "_"
                    + tile_name
                    + aggregate_mode
                    + "_microaggregate_plot_cropped_"
                    + str(image_crop_origin[0])
                    + "_"
                    + str(image_crop_origin[1])
                    + "_"
                    + str(image_crop_size[0])
                    + "x"
                    + str(image_crop_size[1])
                    + "."
                    + fig_format,
                    dpi=300,
                    format=fig_format,
                )
            else:
                plt.savefig(
                    result_dir
                    + "/"
                    + sample_id
                    + "_"
                    + tile_name
                    + aggregate_mode
                    + "tumor_PDL1_immune_PD1_microaggregate_plot"
                    + "."
                    + fig_format,
                    dpi=300,
                    format=fig_format,
                )
            plt.close()

            # save microaggregate size
            for a in microaggregates_final:
                aggregate_sizes.append([sample_id, tile_name, len(a)])

            # aggregate_sizes.append([sample_id, tile_name, len(c)] for c in microaggregates_final)

            # save Tumor PDL1+ / Immune PD1+ counts
            if strict:
                ic_count = coordinates.loc[
                    (coordinates.immune_pheno == 1) & (coordinates.PDL1_marker == 0) & (coordinates.PD1_marker == 1),
                    "immune_pheno",
                ].count()
                tc_count = coordinates.loc[
                    (coordinates.tumor_pheno == 1) & (coordinates.PDL1_marker == 1) & (coordinates.PD1_marker == 0),
                    "tumor_pheno",
                ].count()
            else:
                ic_count = coordinates.loc[
                    (coordinates.immune_pheno == 1) & (coordinates.PD1_marker == 1), "immune_pheno"
                ].count()
                tc_count = coordinates.loc[
                    (coordinates.tumor_pheno == 1) & (coordinates.PDL1_marker == 1), "tumor_pheno"
                ].count()

            cell_counts.append([sample_id, tile_name, ic_count, tc_count])

            #### done with ROI

        # count number of microaggregates in sample
        sample_microaggregate_count = np.sum([len(i) for i in sample_microaggregate_table])

        # rnd_microaggregate_counts = np.asarray(sample_rnd_microaggregate_count_table)
        rnd_microaggregate_counts = np.sum(np.array(sample_rnd_microaggregate_count_table, dtype=object), axis=0)
        # rnd_microaggregate_counts = np.nansum(sample_rnd_microaggregate_count_table)

        # count number of tumor PDL1+ / immune PD1+ paris in sample
        sample_contacts_oi_count = np.sum(sample_contacts_oi_count_table)

        rnd_contacts_oi_counts = np.asarray(sample_rnd_contacts_oi_count_table)

        # Calculate z-score for sample microaggregates
        if np.nansum(rnd_microaggregate_counts) > 0:
            z, p, pr, pl = calculate_zscore_stat(
                sample_microaggregate_count, rnd_microaggregate_counts, num_permutations
            )
        else:
            z = 130
            p = pr = pl = 0.0

        aggregate_counts.append([sample_id, sample_microaggregate_count, z, p, pr, pl])

        # Calculate z-score for sample tumor PDL1+ / immune PD1+ pairs
        if np.nansum(rnd_contacts_oi_counts) > 0:
            z, p, pr, pl = calculate_zscore_stat(sample_contacts_oi_count, rnd_contacts_oi_counts, num_permutations)
        else:
            z = 130
            p = pr = pl = 0.0

        pair_oi_counts.append([sample_id, sample_contacts_oi_count, z, p, pr, pl])

        # Done with sample increase image/sample counter
        i_counter += 1

    # print(aggregate_sizes)
    # print("-----------------")
    # print(aggregate_counts)
    # print("-----------------")
    # print(cell_counts)

    # sys.exit(0)

    df = pd.DataFrame(aggregate_sizes, columns=["sample_id", "ROI", "cluster_size"])
    df.to_csv(result_dir + "/microaggregate_sizes.csv", index=False)

    df = pd.DataFrame(
        aggregate_counts,
        columns=[
            "sample_id",
            "num_microaggregates",
            "z-score",
            "p-value",
            "p-value_right",
            "p-value_left",
        ],
    )
    df.to_csv(result_dir + "/microaggregate_counts.csv", index=False)

    df = pd.DataFrame(
        pair_oi_counts,
        columns=[
            "sample_id",
            "num_pairs",
            "z-score",
            "p-value",
            "p-value_right",
            "p-value_left",
        ],
    )
    df.to_csv(result_dir + "/tumor_PDL1_iummune_PD1_counts.csv", index=False)

    df = pd.DataFrame(cell_counts, columns=["sample_id", "ROI", "ic_count", "tc_count"])
    df.to_csv(result_dir + "/ictc_counts.csv", index=False)
