# MALVIRUS-data
SARS-CoV-2 variant catalog accompanying MALVIRUS.

`get_seqs.sh` updates the NCBI-based SARS-CoV-2 catalog integrated in [MALVIRUS](https://github.com/algolab/malvirus).
Only users that have write access to the [MALVIRUS-data](https://github.com/AlgoLab/MALVIRUS-data) repo (this one) should run this script.

## Instructions

Clone this repository using the `--recursive` flag.

```bash
git clone --recursive git@github.com:AlgoLab/MALVIRUS-data.git
```

Create and activate the conda environment shipped with this repository.

```bash
conda env create -f environment.yaml
conda activate malvirus-data
```

Finally, run `get_seqs.sh`

```bash
./get_seqs.sh
```

The script will compile the required software, download the data from NCBI, build the updated catalog, and create a git commit.

**Please, review the commit before pushing it.**
