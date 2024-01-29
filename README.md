# paiof (& hamtree)

**paiof** - pairwise average amino acid identity using OrthoFinder results (pronounced pay off - hopefully the code isn't corrupt).

**hamtree** - hamming-distance based construction of neighbor-joining tree from orthogroup presence/absence.

![image](https://github.com/raufs/paiof/assets/4260723/cea50cc5-d368-4361-8a6f-09fdc975861c)

Note, these programs work off of OrthoFinder v2.5.4

## Installation:

```shell
# 1. clone Git repo and change directories into it!
git clone https://github.com/raufs/paiof/
cd paiof/

# 2. create conda environment using yaml file and activate it!
conda env create -f paiof_env.yml -n paiof_env
conda activate paiof_env

# 3. complete python installation with the following command:
pip install -e .
```

## Usage

### paoif

```

        Program: paiof
        Author: Rauf Salamzade
        Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

        paiof: Pairwise Average Amino-Acid Identity using OrthoFinder Results.

        Note, currently, the coarse orthogroups are used, not phylogenetic/hierarchical orthogroups.


options:
  -h, --help            show this help message and exit
  -i ORTHOFINDER_RESULTS_DIR, --orthofinder_results_dir ORTHOFINDER_RESULTS_DIR
                        Directory with results from OrthoFinder (Note, this usually starts with Results_DATE/).
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output directory.
  -c CPUS, --cpus CPUS  Number of CPUs to use [Default is 1].
  -scc, --single_copy_core
                        Only use single-copy-core orthogroups.
  -v, --version         Report version of skDER.
```

### hamtree

```
usage: hamtree [-h] -i ORTHOFINDER_RESULTS_DIR -o OUTPUT_DIRECTORY [-p] [-v]

        Program: hamtree
        Author: Rauf Salamzade
        Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology



options:
  -h, --help            show this help message and exit
  -i ORTHOFINDER_RESULTS_DIR, --orthofinder_results_dir ORTHOFINDER_RESULTS_DIR
                        Directory with results from OrthoFinder (Note, this usually starts with Results_DATE/).
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output directory.
  -p, --perform_pca     Also perform PCA analysis based on orthogroup presence/absence.
  -v, --version         Report version of skDER.
```

### Citation

If you find paiof useful for your research, please cite:

* OrthoFinder: phylogenetic orthology inference for comparative genomics. Emms and Kelly, _Genome Biology_ 2016.
* FINAL PAPER NAME, Salamzade, Kalan, and Currie, _bioRxiv_, 2024
