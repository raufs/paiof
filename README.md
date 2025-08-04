# paiof (& hamtree)

**paiof** - pairwise average amino acid identity using OrthoFinder results (pronounced pay-off).

**hamtree** - hamming-distance based construction of neighbor-joining tree from orthogroup presence/absence.

![image](https://github.com/raufs/paiof/assets/4260723/cea50cc5-d368-4361-8a6f-09fdc975861c)

## Dependencies

Note, these programs work off of OrthoFinder v2.5.4/v2.5.5 results. 

Dependencies of paiof & hamtree include:

  - python (tested with v3.10)
  - setuptools
  - pip
  - bioconda::diamond
  - scikit-bio
  - scikit-learn
  - scipy
  - plotly

Specifics on conda environments and versions of additional dependencies used for our manuscript can be found in the subdirectory: `manuscript_analysis_details/`.

## Installation:

Should take < 10 minutes (requires conda installation). Some conda channels with dependencies for paiof only support macOS and Linux.

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

## paiof

### Algorithm Description:

paiof calculates pairwise Average Amino Acid Identity (AAI) between genomes using OrthoFinder ortholog group results. The algorithm works in three main steps:

1. **Orthogroup Processing**: Reads OrthoFinder results to identify orthogroups (either all orthogroups or single-copy core orthogroups if `-scc` flag is used).

2. **DIAMOND Database Creation**: Creates individual DIAMOND databases for each orthogroup's protein sequences to enable efficient similarity searches.

3. **Self-BLAST and AAI Calculation**: 
   - Performs DIAMOND BLASTp searches within each orthogroup (self-blast)
   - Identifies reciprocal best hits (RBH) between genomes within each orthogroup
   - Calculates pairwise AAI as the mean percent identity of all shared orthologous genes between genome pairs
   - Computes standard deviation of AAI values and gene sharing statistics

The algorithm uses checkpoint files to enable resuming interrupted runs and supports multi-threading for parallel processing.

### Example Command:

```
paiof -i orthofinder_results/Results_Aug03/ -o expected_paiof_results/ -c 10
```

The above command will use 10 threads!

### Info on Inputs:

* **OrthoFinder_Results/Results_<MonthDay>/**: Path to the directory of results from OrthoFinder. The tool expects the following subdirectories and files:
  - `Orthogroup_Sequences/`: Directory containing FASTA files for each orthogroup
  - `WorkingDirectory/SpeciesIDs.txt`: Mapping of species IDs to genome names
  - `WorkingDirectory/SequenceIDs.txt`: Mapping of sequence IDs to species IDs
  - `Orthogroups/Orthogroups.tsv`: Orthogroup assignments for all genes
  - `Orthogroups/Orthogroups_SingleCopyOrthologues.txt`: Single-copy orthologue assignments (used with `-scc` flag)

### Info on Outputs:

* **Pairwise_AAI_Estimates.tsv**: Main output file containing pairwise AAI values between all genome pairs with columns:
  - `Genome_1`, `Genome_2`: Genome identifiers
  - `AAI`: Average amino acid identity (mean percent identity of shared orthologs)
  - `Stdev_AAI`: Standard deviation of AAI values
  - `Genes_Shared`: Number of shared orthologous genes
  - `Genes_Shared_Normalized_by_Genome_1_Genes`: Proportion of Genome_1's genes that are shared with Genome_2

* **Full_RBH_Listings.tsv**: Detailed reciprocal best hit information for each gene pair with columns:
  - `Genome_1`, `Genome_2`: Genome identifiers
  - `Gene_1`, `Gene_2`: Gene identifiers
  - `Percent_Identity`: Percent identity of the gene pair

### Usage

```
usage: paiof [-h] -i ORTHOFINDER_RESULTS_DIR -o OUTPUT_DIRECTORY [-c CPUS] [-scc] [-v]

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
  -v, --version         Report version of paoif/hamtree.
```

## hamtree

### Algorithm Description:

hamtree constructs neighbor-joining phylogenetic trees based on orthogroup presence/absence patterns using Hamming distances. The algorithm works as follows:

1. **Orthogroup Presence/Absence Matrix**: Reads the OrthoFinder `Orthogroups.tsv` file to create a binary matrix where each row represents a genome and each column represents an orthogroup. A value of 1 indicates presence of the orthogroup in that genome, while 0 indicates absence.

2. **Hamming Distance Calculation**: Computes pairwise Hamming distances between all genome pairs. Hamming distance measures the proportion of orthogroups that differ in presence/absence between two genomes, providing a measure of genomic content similarity.

3. **Neighbor-Joining Tree Construction**: Uses scikit-bio's neighbor-joining algorithm to construct a phylogenetic tree from the Hamming distance matrix. This tree represents the evolutionary relationships based on shared gene content rather than sequence similarity.

4. **Optional PCA Analysis**: If the `-p` flag is used, performs Principal Component Analysis on the orthogroup presence/absence matrix to visualize genomic content relationships in reduced dimensions.

### Example Command:

```
hamtree -i orthofinder_results/Results_Aug03/ -o expected_hamtree_results/ -p
```

The above command will also perform PCA analysis!

### Info on Inputs:

* **OrthoFinder_Results/Results_<MonthDay>/**: Path to the directory of results from OrthoFinder. The tool expects:
  - `Orthogroups/Orthogroups.tsv`: Orthogroup assignments for all genes, where each row contains an orthogroup ID followed by tab-separated gene lists for each genome

### Info on Outputs:

* **HamTree.newick**: Newick format phylogenetic tree constructed from Hamming distances based on orthogroup presence/absence patterns

* **Plotly_PCA.html** (if `-p` flag used): Interactive HTML file containing PCA visualization of orthogroup presence/absence patterns with:
  - Scatter plot matrix showing relationships between principal components
  - Total explained variance percentage
  - Interactive visualization of genomic content relationships

* **Command_Issued.txt**: Record of the command used
* **Progress.log**: Detailed execution log

### Usage

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
  -v, --version         Report version of paoif/hamtree.
```

### Citation

If you find paiof useful for your research, please cite:

* OrthoFinder: phylogenetic orthology inference for comparative genomics. Emms and Kelly, _Genome Biology_ 2016.
* Complex multicellularity linked with expanded chemical arsenals in microbes. Salamzade, Kalan, and Currie, _in review_, 2025.
