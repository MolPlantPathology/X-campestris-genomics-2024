# Figure 5

## Code run on compute cluster

### Defensefinder

Running `Defensefinder` in a simple for loop:

```bash
(base) [mpaauw@omics-h0 20231206_defensefinder]$ cat defensefinder_runner.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < barcodes.txt

for sample in "${SAMPLES[@]}"
do
	defense-finder run ../20210518_annotation/${sample}/${sample}.faa	
done
```

Defensefinder results were concatenated into one results file `defensefinder_data.tsv`

### CRISPRCasFinder

Run `CRISPRCasFinder` on all genomes to get the spacer arrays. Note that we use `jq` library to filter the resulting `.json` files for only evidence level 4 CRISPR arrays.

```bash
(base) [mpaauw@omics-h0 CRISPRCasFinder]$ cat CRISPR_runner.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < barcodes.txt
awk 'BEGIN {print "barcode", "CRISPRspacers"}' results.tsv

for sample in "${SAMPLES[@]}"
do
	perl CRISPRCasFinder.pl -q -cpuM 12 -in ../../20210317_assembly96/rotated/${sample}_rotated.fasta -out ${sample}_crispr
	
	mv ${sample}_crispr/result.json ${sample}_crispr/${sample}_result.json

    # here we use jq to filter the resulting JSON file for only evidency level 4 CRISPR arrays
	cat ${sample}_crispr/${sample}_result.json | jq 'if has("Crisprs")
    							then .Crisprs |= map(select(.Evidence_Level == 4))
    							else .
    							end
    							| .Date as $date
   							| .Version as $version
    							| .Command as $command
    							| .Sequences |= map(if has("Crisprs") 
								then .Crisprs |= map(select(.Evidence_Level == 4))
                        					else .
                        					end)
    							| .Empty as $empty
    							| {Date: $date, Version: $version, Command: $command, Sequences: .Sequences}' > ${sample}_crispr/${sample}_result_filtered.json
	
	number_of_spacers=$(grep Spacer ${sample}_crispr/${sample}_result_filtered.json | grep -v Spacers | wc -l) 	

	echo -e "$sample\t$number_of_spacers" >> results.tsv
done
```

This is nice, but gives us spacer arrrays of all genomes individually, in a complex `.json` format that is not easy to concatenate into one file for all genomes. It's easier to concatenate all CRISPR containing genomes into one big `fasta` file, and then run the CRISPRCasFinder tool on this concatenated genome file In that way, we get all Xcr CRISPR spacers into one results file!

```bash
(base) [mpaauw@omics-h0 pan_CRISPR]$ cat genome_concatenator.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < CRISPR_plus.txt

touch CRISPR_plus_genomes.fastas

for sample in "${SAMPLES[@]}"
do
	cat genomes/${sample}.fasta >> CRISPR_plus_genomes.fasta
done

# then run CRISPRCasFinder on this file
```

We can paste the pangenome `.json` in [CRISPRTarget](http://crispr.otago.ac.nz/CRISPRTarget/crispr_analysis.html) to find targets. 

### Conservation of each spacer

```bash
TODO
```

In the `R` script to analyse the spacers and targets (see below) we get a list of plasmid IDs that are targetted by the Xcr spacers. We download them using NCBIs efetch

```bash
(base) [mpaauw@omics-h0 20230811_target_plasmids]$ wc -l plasmid_ids.txt 
191 plasmid_ids.txt
(base) [mpaauw@omics-h0 20230811_target_plasmids]$ head -n 4 plasmid_ids.txt 
NZ_CP066960.1
NZ_CP066927.1
NZ_CP066957.1
NZ_LN811400.1

(base) [mpaauw@omics-h0 20230811_target_plasmids]$ cat efetch_batch.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < plasmid_ids.txt

for sample in "${SAMPLES[@]}"
do
	efetch -db nucleotide -format fasta -id ${sample}  > ${sample}.fa
done
```

These plasmid sequences are then annotated using `prokka`, ORF clustered using `roary`, and inspected for ISs using `ISESscan` as before on the full genome set. From the roary results, a selection of T3Es and plasmid related proteins were extracted manually and saved to `xop_selection_presence_absence.csv`.

## Code to generate figures

### Figure 5H and Figure S6
Plot the presence/absence of each defensesystem across all genomes (Figure S6) and the number of defensesystems per genome (Fig 5H) are generated using `defensefinder_plots.R`.

### Figure 5B-F
Plot these panels using `CRISPR_spacers_plots.R`

### Figure S5
Plot this figure using `targetted_plasmid_annotation_analysis.R