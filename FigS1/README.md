# Fig S1

## Calculation of genome characteristics per genome

Here we calculate genome characteristics per genome:

* GC content (note: this was done by actually counting all `C` and `G` bases using `grep`. There are faster ways to do this!)
* Total genome size
* Number of CDS in prokka annotation
* Number of tRNAs in prokka annotation

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

## BUSCO genome completeness

An `sbatch` script to run `BUSCO` on all assemblies after the homopolishing step. This script was manually modified and run 3 times (raw flye assemblies, medaka+racon polished assemblies, homopolished assemblies.)

```bash
(base) [mpaauw@omics-h0 20210322_busco96]$ cat busco_homopolish.sh 
#!/bin/bash
#SBATCH --job-name=busco_96
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=10000
#SBATCH --array=0-95
source activate busco
declare -a SAMPLES
mapfile -t SAMPLES < /zfs/omics/personal/mpaauw/Xc_genomics/00_metadata/barcodes.txt

INPUT="/home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/homopolish/"

echo  "${INPUT}${SAMPLES[$SLURM_ARRAY_TASK_ID]}/consensus_homopolished.fasta" 
busco -f -i "${INPUT}${SAMPLES[$SLURM_ARRAY_TASK_ID]}/consensus_homopolished.fasta" -l xanthomonadales_odb10 -m genome -o "homopolish-${SAMPLES[$SLURM_ARRAY_TASK_ID]}" -c 2 --out_path "homopolish/" --offline
```

## Core genome phylogeny, including *X. campestris* pv. *incanae* genomes

As Figure 1, but includes 5 Xci genomes.