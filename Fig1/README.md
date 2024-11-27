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


Here we calculate genome characteristics per genome:

* GC content
* total genome size
* number of CDS in prokka annotation
* number of tRNAs in prokka annotation

These datasets are plotted in Figure S1.

```bash
(base) [mpaauw@omics-h0 20210819_genome_parameters]$ cat GC_genome_size.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < barcodes.txt

awk 'BEGIN {print "barcode", "genome_size", "GC_content"}' | column -t

for sample in "${SAMPLES[@]}"
do
	total=$(grep -E 'T|A|C|G' -o ../20210317_assembly96/rotated/${sample}_rotated.fasta | wc -l)
	GC=$(grep -E 'C|G' -o ../20210317_assembly96/rotated/${sample}_rotated.fasta | wc -l)	
	echo -n $sample; echo -n -e '\t'; awk "BEGIN {print $total, $GC/$total}" | column -t
done

(base) [mpaauw@omics-h0 20210819_genome_parameters]$ cat prokka_stats.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < barcodes.txt
awk 'BEGIN {print "barcode", "CDS", "tRNA"}' | column -t

for sample in "${SAMPLES[@]}"
do
	CDS=$(grep 'CDS' ../20210518_annotation/${sample}/${sample}.txt | cut -d ' ' -f 2)
	tRNA=$(grep 'tRNA' ../20210518_annotation/${sample}/${sample}.txt | cut -d ' ' -f 2)
	echo -n $sample; echo -n -e '\t'; awk "BEGIN {print $CDS, $tRNA}" | column -t
done
```

## Code to generate figures


