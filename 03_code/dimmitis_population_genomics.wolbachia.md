# Dirofilaria immitis population genomics: genome improvement

### stephen doyle




```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/WOLBACHIA

```


## Porechop to process the raw nanopore data
```bash

bsub.py --threads 10 10 porechop "porechop --input ../GENOME_IMPROVEMENT/NANOPORE_DATA/SRR14299255_1.fastq.gz --output SRR14299255_1_pchoptrim.fastq.gz --threads 10"

```




## Flye assembly of the polished data
```bash
bsub.py --threads 20 50 flye "flye --nano-raw SRR14299255_1_pchoptrim.fastq.gz --genome-size 100M --threads 20 --out-dir flye_out"

# from flye output
#	Total length:	915677
#	Fragments:	1
#	Fragments N50:	915677
#	Largest frg:	915677
#	Scaffolds:	0
#	Mean coverage:	461
```
- single contig with 461x coverage



## Check to see if the assembly maps to the know wolbachia genome 


```bash
minimap2 -x asm20 Dirofilaria_immitis_wolbachia_2.2.fna assembly.fasta


# paf output
#contig_1	915677	7	762601	-	wDi22.scaf1	919954	3	766886	746238	766904	60	tp:A:P	cm:i:127757	s1:i:746081	s2:i:996	dv:f:0.0043	rl:i:0
#contig_1	915677	763413	915677	-	wDi22.scaf1	919954	766893	919947	149308	153054	60	tp:A:P	cm:i:25644	s1:i:149301	s2:i:44	dv:f:0.0043	rl:i:0
```
- clear good match to Wolbachia




## Pilon polish of the Wb genome
- it has been poilished with flye, but will also polish with pilon to be consistent with the original paper 

```bash
module load pilon/1.23-c1

bwa index assembly.fasta

bsub.py --threads 10 10 map "bwa mem -t 10 -x ont2d assembly.fasta SRR14299255_1_pchoptrim.fastq.gz \| samtools sort -o bwa_mapping.sorted.bam"

samtools index bwa_mapping.sorted.bam

bsub.py --threads 10 10 pilon "java -jar ~sd21/lustre_link/software/GENOME_IMPROVEMENT/pilon-1.22/pilon-1.22.jar --genome assembly.fasta --bam bwa_mapping.sorted.bam --threads 10"
```




## mitochondrial genome
- checked for the mitochondrial genome in the ICBAS_JMDir_1.0 however it was missing
- there are no nanopore reads, so will just use one of the two mtDNA genomes from the Di2.2 assembly

```bash

wget http://nematodes.org/downloads/959nematodegenomes/blast/db2/Dirofilaria_immitis_mDi_Athens_and_mDi_Pavia_2.1.fna.gz

samtools faidx Dirofilaria_immitis_mDi_Athens_and_mDi_Pavia_2.1.fna mDi_Athens_2.1 > mDi_Athens_2.1.fa

```


```R
library(ggplot2)

```
