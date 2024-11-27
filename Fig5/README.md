# Figure 5

## Code run on compute cluster

Run CRISPRCasFinder on all genomes to get the spacer arrays. Note that we use `jq` library to filter the resulting `.json` files for only evidence level 4 CRISPR arrays.

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

This is nice, but leaves us with spacer arrrays by all genomes individually, in a complex `.json` format that is not easy to concatenate into one file for all genomes. It's easier to concatenate all CRISPR containing genomes into one big `fasta` file, and then run the CRISPRCasFinder tool on this concatenated genome file In that way, we get all Xcr CRISPR spacers into one results file!

```bash
(base) [mpaauw@omics-h0 pan_CRISPR]$ cat genome_concatenator.sh 
#bin/bash!
while IFS= read -r line; do SAMPLES+=("$line");done < CRISPR_plus.txt

touch CRISPR_plus_genomes.fastas

for sample in "${SAMPLES[@]}"
do
	cat genomes/${sample}.fasta >> CRISPR_plus_genomes.fasta
done

# then run CRISPRTarget on this file
```

We can paste the pangenome `.json` in [CRISPRTarget](http://crispr.otago.ac.nz/CRISPRTarget/crispr_analysis.html) to find targets.


## to be continued

```bash
## add it!
```


## Code to generate figures