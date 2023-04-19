# Dirofilaria immitis populaiton genomics: 
## mapping resistance-associated SNPs

### stephen doyle


## Collecting genomes and SNPs to transfer coordinates
```bash
# working directory

# reference genomes
ln -s ../dimmitis_WSI_2.2.fa .
ln -s ../REFERENCE_GENOMES/dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa .

```

- "old_snp_positions.txt" containing scaffold and position of 42 SNPs
```bash
nDi.2.2.scaf00001	466197
nDi.2.2.scaf00001	748216
nDi.2.2.scaf00004	79159
nDi.2.2.scaf00004	79766
nDi.2.2.scaf00005	662854
nDi.2.2.scaf00007	300005
nDi.2.2.scaf00007	375510
nDi.2.2.scaf00019	442063
nDi.2.2.scaf00021	25243
nDi.2.2.scaf00021	212599
nDi.2.2.scaf00021	379745
nDi.2.2.scaf00021	387398
nDi.2.2.scaf00023	336087
nDi.2.2.scaf00046	22857
nDi.2.2.scaf00046	76278
nDi.2.2.scaf00046	222254
nDi.2.2.scaf00056	97452
nDi.2.2.scaf00056	130632
nDi.2.2.scaf00107	106233
nDi.2.2.scaf00140	30919
nDi.2.2.scaf00185	10639
nDi.2.2.scaf00185	62174
nDi.2.2.scaf00215	28260
nDi.2.2.scaf00238	10209
nDi.2.2.scaf00238	29165
nDi.2.2.scaf00284	33605
nDi.2.2.scaf00284	42920
nDi.2.2.scaf00293	46955
nDi.2.2.scaf00377	15477
nDi.2.2.scaf00492	12704
nDi.2.2.scaf00495	19924
nDi.2.2.scaf00582	14587
nDi.2.2.scaf00589	15334
nDi.2.2.scaf00597	12915
nDi.2.2.scaf00664	25005
nDi.2.2.scaf00669	3266
nDi.2.2.scaf00706	13761
nDi.2.2.scaf01340	1522
nDi.2.2.scaf01422	4176
nDi.2.2.scaf01527	5968
nDi.2.2.scaf06378	423
nDi.2.2.scaf06614	182

```

## Extracting postions from old genome and identifing position in new genome
```bash
# make a bed file of SNPs postions, with an additional coordinate 100 bp upstream 
awk '{print $1,$2,$2+100}' OFS="\t" old_snp_positions.txt > old_snp_positions.100bp.bed
awk '{print $1,$2,$2+900}' OFS="\t" old_snp_positions.txt > old_snp_positions.900bp.bed

# extract fasta sequence of the bed
bedtools getfasta -fi dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa -fo old_snp_positions.100bp.fasta -bed old_snp_positions.100bp.bed
bedtools getfasta -fi dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa -fo old_snp_positions.1000bp.fasta -bed old_snp_positions.1000bp.bed

# use nucmer to map fasta sequences to new genome
nucmer old_snp_positions.100bp.fasta dimmitis_WSI_2.2.fa

show-coords -lTH out.delta | wc -l
#> 40
#> 2 SNPs not found

show-coords -lTH out.delta

```
- transferred coords into the Supplementary data table
- note that the coords were out of sync by one base, likely due to bedtools coord system - sometimes plus sometimes minus.
- confirmed by finding SNPs in the vcffile from Javi "/nfs/users/nfs_s/sd21/team333_link/jg34/POP_Diro/04_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.8.recode.vcf" 


- checking for missing SNPs by mapping 10kb rather than 100 bp
```bash
awk '{print $1,$2,$2+10000}' OFS="\t" old_snp_positions.txt > old_snp_positions.10000bp.bed
bedtools getfasta -fi dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa -fo old_snp_positions.10000bp.fasta -bed old_snp_positions.10000bp.bed
nucmer -p map_10000bp_old_to_new old_snp_positions.10000bp.fasta dimmitis_WSI_2.2.fa
show-coords -lTH map_10000bp_old_to_new.delta
```
- there were no hits to 10kb fragments, which suggests they are truely missing from the genome
- they are big scaffolds in the original assembly, and some parts do map, just not these. 




```bash
samtools faidx dimmitis_WSI_2.2.fa
cat dimmitis_WSI_2.2.fa.fai | grep "chr" | grep -v "chrWb\|chrMtDNA" | sort  > chr.data


for i in `ls *list`; do
vcftools --vcf nuclear_samples3x_missing0.8.recode.vcf --keep ${i} --positions new_snp_postions.txt --extract-FORMAT-info AD --out ${i};
done


for i in ` ls *list.AD.FORMAT`; do 
    cat ${i} | cut -f1,2 | grep -v "CHROM" > pos
    cat ${i} | grep -v "CHROM" | sed 's/,/\t/g' | awk '{sum=0; for (i=3; i<=NF; i+=2) { sum+= $i } print sum}' > sum_ref
    cat ${i}  | grep -v "CHROM" | sed 's/,/\t/g' | awk '{sum=0; for (i=4; i<=NF; i+=2) { sum+= $i } print sum}' > sum_var
    paste pos sum_ref sum_var | awk -v name=${i%.list.AD.FORMAT} '{print name,$0,$4/($3+$4)}' OFS="\t" ; 
    done > snp_freq.pop.data.txt
```

```R
library(tidyverse)
library(reshape2)
library(patchwork)
library(viridis)


chr_data1 <- read.table("chr.data", header=F)
chr_data2 <- chr_data1 %>% mutate(V2=1)
chr_data <- rbind(chr_data1, chr_data2)
chr_data <- chr_data %>% select(V1,V2) %>% mutate(variable="Bourguinat_2017", value="NA") %>% rename(chromosome = V1, position = V2)

data1 <- read.table("snp_study.data.txt",header=T, sep="\t")

data1 <- data1 %>% melt(id.vars=c("chromosome","position"))
data1 <- data1 %>% filter(!is.na(value))

data1 <- rbind(data1, chr_data)

plot_1 <- ggplot() + 
    geom_point(data=data1,aes(position/1e6,variable, alpha=value), size=1, show.legend = FALSE, col="blue") +
    facet_grid(~chromosome, scales="free_x") + 
    theme_bw() +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),
    axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(title="e", x="", y="Study", fill="Allele freq.") +
    scale_alpha_manual(values=c(1,0)) 


data2 <- read.table("snp_freq.pop.data.txt", header=F)

plot_2 <- ggplot(data2) + 
    geom_tile(aes(x=factor(V3), y=V1, fill=V6)) + 
    facet_grid(~V2, scales="free_x", space="free_x") +
    scale_fill_gradient(low="white", high="blue", limits = c(0, 1)) +
    #scale_fill_viridis(option="magma", direction=-1) +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="bottom") +
    labs(x="SNP position", y="Country", fill="Allele freq.")


plot_1 + plot_2 + plot_layout(ncol=1)
```



## Rerunning Javi's Dxy data to make a plot that integrates with the resistance panels
- Javi shared the file "country_dxy.txt"

```R
library(tidyverse)
library(reshape2)
library(patchwork)
library(viridis)


data3 <- read.table("country_dxy.txt", header=T)
data3 <- data3 %>% mutate(comparison = paste0(pop1,"_vs_",pop2))

plot_3 <- ggplot(data3) +
    geom_line(aes(x=((window_pos_1+window_pos_2)/2)/1e6, y=avg_dxy, col=comparison), size=.5) + facet_grid(~chromosome, scales="free_x") + theme_bw() + 
    scale_color_brewer(type = 'div', palette = 'Accent') +
    labs(title="d", y="Dxy", x="Genomic position (Mb)") +
    theme(legend.position="top")



# incorporating plots from above
plot_3 + plot_1 + plot_2 + plot_layout(ncol=1, height=c(2,1,1))

ggsave("Dxy_resistant_snps.genomewide.pdf", height=7, width=10)
ggsave("Dxy_resistant_snps.genomewide.png")
```