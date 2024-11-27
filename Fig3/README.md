## Code run on compute cluster

### ISEScan and distance between ISs and T3Es
`sbatch` script to run ISESscan on all genomes:

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

Then we calculate the distance between each T3E to the closest IS in each genome:

```bash
TODO
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

### Structural variation caller

```bash
TODO (this is a lot of code)
```

