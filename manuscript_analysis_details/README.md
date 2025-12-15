## Reproducing Results from paiof and hamtree in the Associated Manuscript

To reproduce results from hamtree's usage in the manuscript, please first run OrthoFinder (v2.5.4) on the input set of 133 Actinomycetota proteomes (after gunzip-ing them). Details on the OrthoFinder conda environment can be found in this current directory: `conda_env_with_orthofinder.yaml`. 

Afterwards, you would just run hamtree (v1.0.0) using an environment similar/identical to the one described in `paiof_conda_env.yaml`:

```
hamtree -i OrthoFinder_Results/Results_<MonthDay>/ -o Hamtree_Results/ -p
```

Note, in the revised manuscript, we only use paiof for a supplemental analysis that is briefly described in the supplemental methods document. The details of this analysis, performed using paiof (v1.0.2) are included here, with the final plot generated included as a supplemental figure in the GitHub repo associated with the manuscript at: https://github.com/Kalan-Lab/Ancestral_BGC_Exansions_Study

Info on the environment used to run paiof can be found in the file `paiof_conda_env.yaml`. 
