## Test Case

Here we include a test case for running paiof and hamtree analysis on results from OrthoFinder (v2.5.4) pertaining to 4 samples (proteomes for which can be found in FASTA format in the `proteomes/` subdirectory).

> [!NOTE]
> Expected results for paiof and hamtree can be found in the compressed folders: `expected_paiof_results.tar.gz` and `expected_hamtree_results.tar.gz`, respectively. If you run and compare output files to expected output they might differ in the ordering of lines, but if you use the `sort` command to reorder lines alphabetically you should see that the contents are identical.

### Run OrthoFinder

To create the required input for paiof and hamtree, we must first run OrthoFinder (either v2.5.4 or v2.5.5 - will update future versions to also handle OrthoFinder v3+ results).

```
orthofinder -f input_proteomes/ -o orthofinder_results/ -c 20  
```

> The timing of the above command using 20 threads should be around 5.07 minutes.

> [!NOTE]
> You can skip this step and just use the provided OrthoFinder results for this test case.

### Run paiof

Now we can run paiof to perform intra-ortholog group pairwise alignments and report pairwise AAI between samples.

```
paiof -i orthofinder_results/Results_Aug03/ -o expected_paiof_results/ -c 20
```

> The timing of the above command using 20 threads should be around 14.45 minutes.

### Run hamtree 

Now we can run hamtree to create a neighbor-joining tree constructed from Hamming distances between samples based on ortholog group carriage profiles.

```
hamtree -i orthofinder_results/Results_Aug03/ -o expected_hamtree_results/
```

> The timing of the above command using 20 threads should be around 1.75 minutes.
