## Code run on compute cluster

### ISEScan
`sbatch` script to run ISESscan on all genomes. The number of insertion sequences per genome is plotted across the phylogenetic tree in Figure 3C, and as boxplots divided by Xcc and Xcr in Figure S2. 

```bash
(base) [mpaauw@omics-h0 20211007_transposons]$ cat ISESscan_runner.sh 
#!/bin/bash
#
#SBATCH --job-name=isescan_96
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=6000
#SBATCH --array=0-95

source activate isescan

declare -a SAMPLES
mapfile -t SAMPLES < barcodes.txt

INPUT="/home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/rotated/"

prefix=${SAMPLES[$SLURM_ARRAY_TASK_ID]} 
prefix="${prefix/barcode/bc}"

isescan.py --seqfile "${INPUT}${SAMPLES[$SLURM_ARRAY_TASK_ID]}_rotated.fasta" --output "${SAMPLES[$SLURM_ARRAY_TASK_ID]}" --nthread 8
```

### MobileOGdb

```bash
(base) [mpaauw@omics-h0 20230926_mobileOGdb]$ cat mobileOG_runner.sh 
#!/bin/bash
#
#SBATCH --job-name=isescan_96
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=6000
#SBATCH --array=0-95

source activate mobileOG-db 

declare -a SAMPLES
mapfile -t SAMPLES < barcodes.txt

bash mobileOGs-pl-kyanite.sh -i genomes/${SAMPLES[$SLURM_ARRAY_TASK_ID]}_rotated.fasta -d mobileOG-db-beatrix-1.6.dmnd -m mobileOG-db-beatrix-1.6-All.csv  -k 15 -e 1e-15 -p 80 -q 80
```

### Identifying regions of genome plasticity using ppanggolin

Here, we rerun ppanggolin using only the Xcr/Xcc genomes (so without the outgroup and nonpathogenic isolate) to identify regions of genome plasticity (RGPs) between Xcc and Xcr genomes.

```bash
ppanggolin workflow --anno fasta-list.txt
ppanggolin rgp -p pangenome.h5
ppanggolin write -p pangenome.h5 --regions --output RGPs

(base) [mpaauw@omics-h0 RGPs]$ head -n 4 plastic_regions_v2.tsv 
region	organism	contig	start	stop	genes	contigBorder	wholeContig
Xcc_barcode05_contig_1_polish_RGP_6	Xcc_barcode05	contig_1_polish	55852	83600	23	False	False
Xcc_barcode05_contig_1_polish_RGP_18	Xcc_barcode05	contig_1_polish	161542	169233	10	False	False
Xcc_barcode05_contig_1_polish_RGP_28	Xcc_barcode05	contig_1_polish	280048	289269	7	False	False
```

We used this file to count in each genome the number of RGPs as well as the number of genes contained within those regions of genome plasticity.

### Structural variation caller

First, we use MUM&Co to call structural variants between all pairwise genome comparisons:

```bash
(base) [mpaauw@omics-h0 20230904_SVs]$ head -n 4 combinations.txt 
barcode01 barcode01
barcode01 barcode02
barcode01 barcode03
barcode01 barcode04

(base) [mpaauw@omics-h0 MUMandCo]$ cat MUM_batch.sh 
#!/bin/bash
#SBATCH --job-name=mum
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-215

declare -a SAMPLES
mapfile -t SAMPLES < comb_split_09  # Note that SBATCH arrays can hold maximum 1000 jobs, therefore we had to split the combinations file in 9 batches.

INPUT="/home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/rotated/"

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]} 
reference="${sample:0:9}"
query="${sample:10:9}"
jobname="${reference}_${query}_output"

mkdir jobs/${jobname}
cd jobs/${jobname}
bash ../../mumandco_v3.8.sh -r /home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/renamed_contigs/${reference}_rotated.fasta -q /home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/renamed_contigs/${query}_rotated.fasta -o ${jobname} -g 5000000 -t 1
```

Then, all ~9000 results file are concatenated into one major structural variant database.

Each genome inversion has two breakpoints: the locations on the genome where the inversion starts, and ends. To calculate the number of inversion breakpoints per genome, we
first split row (each inversion event) into two rows (each inversion breakpoint)
 keep track of each inversion breakpoint, and 'merge' them when they occur on the same spot, or within 150 bp. 

```python
"""
This script gets the structural variant inversions tsv
It select the chromosome, and splits each inversion into two 'breakpoint events' (i.e., the start and end of the inversion)
TODO - output needs to be seperated by barcode to get the inversion breakpoints per barcode
post processing: change contig name to 'contig_1_polish' to be compatible with older data
post processing2: use bedtools merge to merge breakpoints that are very close together. 
"""
file_path = 'SV_inversions_duplicates_removed.tsv'  # Replace with the path to your file
with open(file_path, 'r') as file:
	for line in file:
		parts = line.strip().split('\t')  # Use .strip() to remove leading/trailing whitespac
	
		query = parts[0]
		q_start = int(parts[2])
		q_end = int(parts[3])

		print(query, q_start-1, q_start, sep='\t') # breakpoint 1
		print(query, q_end-1, q_end, sep='\t') #breakpoint 2
```

Then we seperate this information by barcode, and merge inversion breakpoints occuring at a distance within 50 bp into one breakpoint.

```bash
#!/bin/bash

# Specify the path to your .txt file
file_path=barcodes.txt

# Loop over the lines in the file
while IFS= read -r sampleID; do
    # Do something with each sampleID (replace with your desired operation)
    echo "Processing sampleID: $sampleID"

    grep $sampleID inversion_breakpoints_database_duplicates_removed.bed | sort -k1,1 -k2,2n | bedtools merge -d 50 > ${sampleID}_inversion_breakpoints.bed    
done < "$file_path"

```

The number of inversions between gall pairwise genome comparisons is shown in Figure S3A, and number of inversion breakpoints in plotted in figure S3B. The overlap between inversion breakpoints and insertion sequences is shown in Figure 3.

```bash
#!/bin/bash
file_path=barcodes.txt

while IFS= read -r sampleID; do
	# Do something with each sampleID
	cat genome_files/${sampleID}.genome | sed 's/_polish//g' > tmp.genome	
	grep insertion ISESscan_gffs/${sampleID}_rotated.fasta.gff | sed 's/_polish//g' > tmp.IS.gff
	cat inversion_breakpoints_persample/${sampleID}_inversion_breakpoints.bed | sed "s/${sampleID}_//g" | bedtools slop -b 50 -g tmp.genome > tmp.inv.bed

	n_IS=$(cat tmp.IS.gff | wc -l)
	n_overlap=$(bedtools intersect -u -a tmp.IS.gff -b tmp.inv.bed | wc -l)

	total=0

	# shuffle 25 times and calculate the mean of the randomized overlaps
	for ((i=1; i<=25; i++))
	do
		bedtools shuffle -i tmp.IS.gff -g tmp.genome > tmp.shuffle
		n_random=$(bedtools intersect -u -a tmp.shuffle -b tmp.inv.bed | wc -l)
		total=$((total+$n_random))	
	done

	mean=`echo  $total / 25 | bc -l`  
	echo ${sampleID} ${n_IS} ${n_overlap} ${mean} 

done < "$file_path"
```


## Code for plots

See `figure_3_scripts.R`, `mobileOG_plot.R`, `IS_breakpoint_overlaps.R`.