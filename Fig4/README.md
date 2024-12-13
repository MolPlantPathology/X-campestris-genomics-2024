# Figure 4

## Code run on compute cluster

### GC content and codon adaptation index per gene

The calculations for this section are all done using command line utilities in the `EMBOSS` suite, and were performed on our inhouse assemblies/annotations of the reference strains Xcc8004 (barcode71) and Xcr756 (barcode01). To calculate the GC content and codon adaptation per gene, we need seperate `.fasta` files for each gene in both genomes. This was done using `seqretsplit`. Then, we calculate genomewide codon usage statistics on the original `.ffn` containing all open reading frames using `cusp`:
 
```bash
cusp -sequence barcode71.ffn -outfile barcode71.cusp
```

Then, for each gene the codon adaptation index is calcuted using `cai` (using the genome wide codon usage statistics as a reference). Similarly, the GC content is calculated for each gene using `geecee`

```bash
# codon adaptation index
for file in *.fasta
do
	cai -seqall $file  -cfile barcode71.cusp -stdout y -auto y
done

# GC content
for file in *.fasta
do
	geecee $file -stdout y -auto y
done
```

Result files were concatenated and plotted in R (see code to generate figures below)

### Distance between effectors and insertion sequences

For Figure 4B, we calculate the distance between each T3E - insertion sequence pair using `bedtools closest`. We do this for the whole pangenome, in other words, the calculations are performed for each genome assembly and saved per genome. After that, the results were concatenated into a total distance file `distance_total.tsv`.

```bash
(emboss) [mpaauw@omics-h0 20211007_transposons]$ cat distance_T3E_IS.sh 
#bin/bash!
counter=$((1))
while IFS= read -r line; do SAMPLES+=("$line");done < barcodes.txt

for sample in "${SAMPLES[@]}"
do
	awk 'NR>2' ${sample}/rotated/${sample}_rotated.fasta.out | column -t | awk '{print $1" "$2 " " $3 " " $4 " " $5}' | while read line
	do
		contig=$(echo $line | cut -f1 -d" ")
		IS=$(echo $line | cut -f2 -d" ")
		ISs=$(echo $line | cut -f3 -d" ")	
		star=$(echo $line | cut -f4 -d" ") 
		end=$(echo $line | cut -f5 -d" ") 

		echo -e $contig'\t'$star'\t'$end'\t'$IS  >> ${sample}_ISs.bed	
	done
	grep -E 'Xop|Avr|Hop|Rip' ../20210518_annotation/${sample}/${sample}.gff | cut -f1,4,5,9 > ${sample}_T3E.gff
	bedtools closest -d -a ${sample}_T3E.gff -b ${sample}_ISs.bed > ${sample}_distance.bed

	sed -i "s/$/\t${sample}/" ${sample}_distance.bed

	mv ${sample}_distance.bed distances/
	
	rm ${sample}_ISs.bed
	rm ${sample}_T3E.gff	
done
```

For Figure 4C, we inspect only our assembly of Xcc8004 (barcode71). In this genome, we performed randomization of T3E location using `bedtools shuffle` and calculate again the distance between each randomly placed T3E and the closest insertion sequence.

```bash
(bedtools) [mpaauw@omics-h0 2024_revision_experiment]$ cat bedtools_closest.sh 
bedtools closest -d -a barcode01_T3E.gff -b barcode01_ISs.bed | cut -f 4,9 | awk 'BEGIN {OFS="\t"} {print "barcode01", "Observed", "Trial1", $0}' | sed 's/ID=.*product=//' > closest_result.tsv

for i in {1..100}
do
   bedtools shuffle -i barcode01_T3E.gff -g barcode01.genome > bc01.shuf.t3e.gff
   sort -k1,1 -k2,2n bc01.shuf.t3e.gff > bc01.shuf.t3e.sorted.gff
   bedtools closest -d  -a bc01.shuf.t3e.sorted.gff -b barcode01_ISs.bed | cut -f 4,9 |  awk -v iter="$i" 'BEGIN {OFS="\t"} {print "barcode01", "Randomized", "Trial" iter, $0}' | sed 's/ID=.*product=//' >> closest_result.tsv
done

bedtools closest -d -a barcode71_T3E.gff -b barcode71_ISs.bed | cut -f 4,9 | awk 'BEGIN {OFS="\t"} {print "barcode71", "Observed", "Trial1", $0}' | sed 's/ID=.*product=//' >> closest_result.tsv

for i in {1..100}
do
   bedtools shuffle -i barcode71_T3E.gff -g barcode71.genome > bc71.shuf.t3e.gff
   sort -k1,1 -k2,2n bc71.shuf.t3e.gff > bc71.shuf.t3e.sorted.gff
   bedtools closest -d  -a bc71.shuf.t3e.sorted.gff -b barcode71_ISs.bed | cut -f 4,9|  awk -v iter="$i" 'BEGIN {OFS="\t"} {print "barcode71", "Randomized", "Trial" iter, $0}' | sed 's/ID=.*product=//' >> closest_result.tsv
done
```

## Code to generate figures

### GC content and codon adaptation index per gene

See `GC_CAI_plots.R` for Figure 4A.

### Distance between effectors and insertion sequences

See `distance_T3E_transposon.R` for Figure 4B. See `distance_T3E_transposon_shuffled.R.R` for Figure 4C.