# Figure 1

## Code run on compute cluster

### Pangenome computation: Roary

Running pangenome software software `Roary`. As input, we take the `prokka` derived `.gff` files. resulting in a gene-presence absence matrix. This matrix is later used to count the number of Type III Effector homologs in each genome by searching for the words `Xop` `Hop` or `Rip`. This was done manually in Excel.

```bash
roary -p 24 -e -n -r -v /home/mpaauw/personal/Xc_genomics/04_results/20210518_annotation/barcode*/barcode*.gff
```
### Core genome alignment and phylogenetic tree: PPanGGOLiN & FastTree

We use a different pangenome tool `PPanGGOLiN` to make a core gene alignment, which is later used by `FastTree` to make a phylogenetic tree. The phylogenetic tree is rooted and visualised later in R.

```bash
ppanggolin workflow --anno fasta-list.txt # with fasta-list.txt a tab seperated file with strain identifiers, and the location of the .gff annotation files
ppanggolin msa -p pangolin.h5 –source dna –phylo

# cd into the directory where the core genome alignment is stored
fasttree -nt -gtr core_genome_alignment.aln > core_gene_tree.newick 
```

To calculate pairwise average nucleotide identity values between all genomes, we used `fastANI`. Here, `assemblies.txt` is simply a file specifying the path to all assemblies.

```bash
fastANI --ql assemblies.txt --rl assemblies.txt -t 32 -o all_vs_all.tsv
```

## Code to generate figures

### Bioassay data and annotated phylogenetic tree

See file `figure_1_scripts_v2.R`

### FastANI plot

See file `Figure_1_fastANI_panel.R`
