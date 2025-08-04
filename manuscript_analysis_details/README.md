## Reproducing Results from paiof and hamtree in the Associated Manuscript

To reproduce results from paiof and hamtree used in the manuscript, please first run OrthoFinder (v2.5.4) on the input set of 133 Actinomycetota proteomes (after gunzip-ing them). Details on the OrthoFinder conda environment can be found in this current directory: `conda_env_with_orthofinder.yaml`. 

Afterwards, you would just run paiof (v1.2.0) using an environment similar/identical to the one described in `paiof_conda_env.yaml`:

```
paiof -i OrthoFinder_Results/Results_<MonthDay>/ -o paiof_Results/ -c 20
```

To reproduce hamtree results for the manuscript, please run the following in the same conda environment:

```
hamtree -i OrthoFinder_Results/Results_<MonthDay>/ -o Hamtree_Results/ -p
```

Select paiof results from the manuscript analysis are included in the subdirectory `paiof_Results/`. Similarly, hamtree results used for the manuscript can be found in the subdirectory: `Hamtree_Results/`. 