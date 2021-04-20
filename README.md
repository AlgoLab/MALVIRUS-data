# MALVIRUS-data
SARS-CoV-2 variant catalog accompanying MALVIRUS.

The workflow `Snakefile` updates the NCBI-based SARS-CoV-2 catalog integrated in [MALVIRUS](https://github.com/algolab/malvirus).
Only users that have write access to the [MALVIRUS-data](https://github.com/AlgoLab/MALVIRUS-data) repo (this one) should run this script.

## Instructions

Clone this repository:

```bash
git clone git@github.com:AlgoLab/MALVIRUS-data.git
```

Create and activate the conda environment shipped with this repository:

```bash
conda env create -f environment.yml
conda activate malvirus-data-snakemake
```

Finally, run Snakemake

```bash
make -C software
snakemake -j4 --use-conda -pr
```

The script will compile the required software, download the data from NCBI, build the updated catalogs.

Then, commit the results.


