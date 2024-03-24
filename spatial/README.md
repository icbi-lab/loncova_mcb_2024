# Spatial analyses of IMC

This repo conatins an updated version of the scripts used in [Plattner et al. *iScience*. 2023](https://doi.org/10.1016/j.isci.2023.108399)

We suggest to create a dedicated conda environment and install all software dependencies to run the IMC analysis
python scripts therein.

You can create a conda environment as follows

```bash
conda env create --name IMCanalysis python=3

conda activate IMCanalysis
pip install pathlib
pip install scipy
pip install tqdm
pip install matplotlib
pip install networkx
pip install pandas
pip install seaborn
```

## Obtain the data
All preprocessed data required to run these analyses are available from [Zenodo](https://doi.org/10.5281/zenodo.7015015).

```bash
cd ../
wget -O imaging.zip 'https://zenodo.org/record/7015015/files/imaging.zip?download=1'
unzip imaging.zip
cd spatial
```

## Run the analyses
To run the IMC voronoi analysis as described in the MCB chapter run the following commands:

```bash

mkdir -p ../results/spatial

./spatialAnalysis_voronoi.py \
    --result_dir ../results/spatial \
    --data_dir ../imaging/Single_cell_data \
    --sample_sheet sample_sheet.tsv \
    --celltype_def phenotype_colors.tsv\
    --num_cores 4 \
    --legend True \
    --skip_phenotype Unknown CD57+_cells CD38+_cells

./plotVoronoi.py \
    --result_dir ../results/spatial \
    --data_dir ../imaging/Single_cell_data \
    --sample_sheet sample_sheet.tsv \
    --celltype_def phenotype_colors.tsv 

./cohort_zscore_heatmap.py  \
    --result_dir ../results/spatial \
    --data_dir ../results/spatial \
    --sample_sheet sample_sheet.tsv \
    --celltype_def phenotype_colors.tsv
    --s_count \
    --direction both \
    --discrete_cm \
    --skip_phenotype Unknown CD57+_cells CD38+_cells

./microaggregates_PD1_PDL1.py \
    --result_dir ../results/spatial \
    --data_dir ../imaging/PD1_PDL1 \
    --sample_sheet sample_sheet.tsv \
    --celltype_def phenotype_colors.tsv \
    --num_cores 4

```


## Contact
If you have questions regarding the analyses, please use the [issue tracker](https://github.com/icbi-lab/loncova_mcb_2024/issues).

