# Figure 1

## Code run on compute cluster

Running pangenome software software `roary`, resulting in a gene-presence absence matrix. This matrix is later used to count the number of Type III Effector homologs in each genome by searching for the words `Xop` `Hop` or `Rip`. This was done manually in Excel.

```bash
roary -p 24 -e -n -r -v /home/mpaauw/personal/Xc_genomics/04_results/20210518_annotation/barcode*/barcode*.gff
```

We use a different pangenome tool `ppanggolin` to make a core gene alignment, which is later used by `FastTree` to make a phylogenetic tree

```bash
## add it!
```

To calculate pairwise average nucleotide identity values between all genomes, we used `fastANI`

```bash
## add it!
```


## Code to generate figures


