# Genome assembly and annotation

We sequenced *X. campestris* genomic DNA using Oxford Nanopore Technologies long read technologies. The raw reads are available in NCBI SRA with BioProject number `PRJNA1104797`. Assembly was done using `Flye`, followed by a sequential polishing pipeline using `Racon`, `Medaka`, and later, `Homopolish`.

## Code run on compute cluster

The first three steps were run using [Snakemake](https://snakemake.github.io) - a workflow tool to chain together a set of commands that depend on each others outcome. For this we use the Snakefile and `cluster.json`configuration file present in this repository.

```bash
snakemake -n -j 4 --latency-wait 60 --use-conda --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -c {cluster.c} -t {cluster.time} --mem {cluster.memory} -o {cluster.output}"
```
Then, `homopolish` was run using a `sbatch` script:

```bash
#!/bin/bash
#BATCH --job-name=homopolish_96
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=10000
#SBATCH --array=1-95

declare -a SAMPLES
mapfile -t SAMPLES < /zfs/omics/personal/mpaauw/Xc_genomics/00_metadata/barcodes.txt

INPUT="/home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/medaka_polishing/${SAMPLES[$SLURM_ARRAY_TASK_ID]}/consensus.fasta"

python3 homopolish.py polish -t 2 -a "${INPUT}" -s bacteria.msh -m R9.4.pkl -o "${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
```

And similary rotated using a `sbatch` script to run `Circlator`

```bash
#!/bin/bash
#BATCH --job-name=circlator_96
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=10000
#SBATCH --array=0-95
source activate circlator

declare -a SAMPLES
mapfile -t SAMPLES < /zfs/omics/personal/mpaauw/Xc_genomics/00_metadata/barcodes.txt

INPUT="/home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/homopolish/${SAMPLES[$SLURM_ARRAY_TASK_ID]}/consensus_homopolished.fasta"

circlator fixstart "$INPUT" "${SAMPLES[$SLURM_ARRAY_TASK_ID]}_rotated"  
```

The resulting genome assemblies are also available on NCBI (`PRJNA1104797`), or available as the 'inhouse' version at [Figshare (University of Amsterdam)](https://doi.org/10.21942/uva.25359736.v1). Then, the genomes were annotated using `Prokka` with a home-made database with known virulence factors of plant pathogens, containing 5432 entries (`db_unique.fa`).

```bash
(base) [mpaauw@omics-h0 20210518_annotation]$ cat prokka_scratch.sh 
#!/bin/bash
#
#SBATCH --job-name=prokka_96
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=10000
#SBATCH --array=0-95
source activate prokka

declare -a SAMPLES
mapfile -t SAMPLES < /zfs/omics/personal/mpaauw/Xc_genomics/00_metadata/barcodes.txt

INPUT="/home/mpaauw/personal/Xc_genomics/04_results/20210317_assembly96/rotated/"

prefix=${SAMPLES[$SLURM_ARRAY_TASK_ID]} 
prefix="${prefix/barcode/bc}"

cd /scratch
cp "/zfs/omics/personal/mpaauw/Xc_genomics/04_results/20210518_annotation/db/db_unique.faa" ./

prokka --force --proteins "db_unique.faa" --locustag "$prefix" --outdir "${SAMPLES[$SLURM_ARRAY_TASK_ID]}" --prefix "${SAMPLES[$SLURM_ARRAY_TASK_ID]}" --cpus 8 "${INPUT}${SAMPLES[$SLURM_ARRAY_TASK_ID]}_rotated.fasta"

cp -r "${SAMPLES[$SLURM_ARRAY_TASK_ID]}/" "/zfs/omics/personal/mpaauw/Xc_genomics/04_results/20210518_annotation/${SAMPLES[$SLURM_ARRAY_TASK_ID]}" 

```

## Data availability

The raw reads and genome assemblies are available in NCBI under BioProject `PRJNA1104797`
