
# Dirofilaria immitis populaiton genomics: 
## Comparative Genomics

### stephen doyle



### Dotplot comparing Di and Bm
```bash
grep ">" brugia_malayi.PRJNA10729.WBPS17.genomic.fa | grep "Chr" | sed 's/>//'g | while read -r NAME; do 
    samtools faidx brugia_malayi.PRJNA10729.WBPS17.genomic.fa ${NAME} >> bm_chr.fa; 
    done

grep ">" dimmitis_WSI_1.0.fa | grep "chr" | sed 's/>//g' | while read -r NAME; do 
    samtools faidx dimmitis_WSI_1.0.fa ${NAME} >> di_chr.fa; 
    done


minimap2 -x asm20 bm_chr.fa di_chr.fa > di_v_bm.paf


cat di_v_bm.paf | /nfs/users/nfs_s/sd21/lustre_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout > di_v_bm.layout
cat di_v_bm.layout | cut -f1-19 > di_v_bm.layout2
```


```R
library(tidyverse)
library(viridis)


data<-read.table("di_v_bm.layout2")
data<-data[data$V17 > 300, ]
tdata<-data[data$V17 > 300,  ]
vdata<-aggregate(data$V5, by=list(data$V4), max)

ggplot()+
     geom_segment(data=data,mapping=aes(y=V10/1e6, yend=V11/1e6, x=V15/1e6, xend=V16/1e6, colour=V4)) +
     theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
     geom_point(data=tdata, aes(x=V15/1e6, y=V10/1e6, colour=V4), size=1)+
     geom_point(data=tdata, aes(x=V16/1e6, y=V11/1e6, colour=V4), size=1)+
     labs(x="Brugia malayi chromosomes (Mb)", y="Dirofilaria immitis chromosomes (Mb)", colour="Chromosome")+
     scale_colour_viridis(discrete = TRUE) + facet_grid(V1~V4, scales="free")

ggsave("di_v_bm.chromosomes.dotplot.pdf")
```









## Circos plot comparing Di and Bm
```bash 
samtools faidx di_chr.fa
samtools faidx bm_chr.fa

cat di_chr.fa.fai bm_chr.fa.fai | awk '{print $1, "1", $2}' OFS="\t" > chromosome_lengths.txt
cat chromosome_lengths.txt > chromosome_lengths.sorted.txt


```
- plot
https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
- links
https://jokergoo.github.io/circlize_book/book/genomic-plotting-region.html#genomic-links
```R
library(tidyverse)
library(circlize)
library(viridis)

chr <- read.table("chromosome_lengths.sorted.txt", header=F)
colnames(chr) <- c("chr", "start", "end")
data<-read.table("di_v_bm.layout2")
data<-data %>% filter(V18>5000)

di_links <- data %>% select(V8, V10, V11, V18)
colnames(di_links) <- c("chr", "start", "end", "value")
bm_links <- data %>% select(V13, V15, V16, V18)
colnames(bm_links) <- c("chr", "start", "end", "value")

colours <- viridis(5)
palette(colours)

grid.col = c(dirofilaria_immitis_chr1 = colours[1], dirofilaria_immitis_chr2 = colours[2], dirofilaria_immitis_chr3 = colours[3], dirofilaria_immitis_chr4 = colours[4], dirofilaria_immitis_chrX = colours[5],
    Bm_v4_Chr3_scaffold_001 = "grey", Bm_v4_Chr1_scaffold_001 = "grey", Bm_v4_ChrX_scaffold_001 = "grey", Bm_v4_Chr4_scaffold_001 = "grey", Bm_v4_Chr2_contig_001 = "grey")


pdf()
circos.clear()

circos.genomicInitialize(chr)
circos.track(ylim = c(0, 1), bg.col = grid.col, bg.border = NA, track.height = 0.05)
circos.genomicLink(di_links, bm_links, border = NA,  col = as.factor(di_links$chr))

circos.clear()
dev.off()
```


### COmparison on new Di chromosomes and ICBAS_JMDir_1.0 assembly
```bash 
minimap2 -x asm20 di_chr.fa ../REFERENCE_GENOMES/GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa > di_chr.v.ICBAS_JMDir_1.0.paf

cat di_chr.v.ICBAS_JMDir_1.0.paf | /nfs/users/nfs_s/sd21/lustre_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout > di_chr.v.ICBAS_JMDir_1.0.layout
cat di_chr.v.ICBAS_JMDir_1.0.layout | cut -f1-19 > di_chr.v.ICBAS_JMDir_1.0.layout2

samtools faidx di_chr.fa
samtools faidx ../REFERENCE_GENOMES/GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa

cat di_chr.fa.fai ../REFERENCE_GENOMES/GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa.fai | awk '{print $1, "1", $2}' OFS="\t" > chromosome_lengths.txt

```







```bash
minimap2 -x asm20 di_chr.fa Onchocerca_volvulus_v3.fa > di_v_ov.paf

cat di_v_ov.paf | /nfs/users/nfs_s/sd21/lustre_link/software/GENOME_ASSEMBLY/minimap/utils/bin/layout > di_v_ov.layout

cat di_v_ov.layout | cut -f1-19 > di_v_ov.layout2

samtools faidx di_chr.fa
samtools faidx Onchocerca_volvulus_v3.fa

cat di_chr.fa.fai Onchocerca_volvulus_v3.fa.fai | awk '{print $1, "1", $2}' OFS="\t" > di_v_ov.chromosome.lengths.txt


# extract only the chromosome parts
grep "OVOC.OM" di_v_ov.chromosome.lengths.txt > di_v_ov.chromosome.lengths.txt2
grep "OVOC.OM" di_v_ov.layout | cut -f1-19 > di_v_ov.layout2
```


```R
library(tidyverse)
library(circlize)
library(viridis)

chr <- read.table("di_v_ov.chromosome.lengths.txt2", header=F)
colnames(chr) <- c("chr", "start", "end")
data<-read.table("di_v_ov.layout2")
data<-data %>% filter(V18>5000)

di_links <- data %>% select(V8, V10, V11, V18)
colnames(di_links) <- c("chr", "start", "end", "value")
ov_links <- data %>% select(V13, V15, V16, V18)
colnames(ov_links) <- c("chr", "start", "end", "value")

colours <- viridis(5)
palette(colours)

grid.col = c(dirofilaria_immitis_chr1 = colours[1], dirofilaria_immitis_chr2 = colours[2], dirofilaria_immitis_chr3 = colours[3], dirofilaria_immitis_chr4 = colours[4], dirofilaria_immitis_chrX = colours[5],
    OVOC.OM1a = "grey", OVOC.OM1b = "grey", OVOC.OM2 = "grey", OVOC.OM3 = "grey", OVOC.OM4 = "grey", OVOC.OM5 = "grey")


pdf("di_v_ov.circlize.pdf")
circos.clear()

circos.genomicInitialize(chr)
circos.track(ylim = c(0, 1), bg.col = grid.col, bg.border = NA, track.height = 0.05)
circos.genomicLink(di_links, ov_links, border = NA,  col = as.factor(di_links$chr))

circos.clear()
dev.off()
```






## Genome wide coverage plots of male and female worms
### Generate coverage data
```bash


```


### plots
```R
library(tidyverse)
library(viridis)
library(patchwork)

#female
data_f <- read.table("/nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/POPGEN/MAPPING_mapping_dimmitis/ERR034940_ERR034940_mapping/ERR034940.10000_window.cov", header=F)

data_m <- read.table("/nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/POPGEN/MAPPING_mapping_dimmitis/SRR13154015_SRR13154015_mapping/SRR13154015.10000_window.cov", header=F)



data_f <- data_f %>% filter(str_detect(V1, 'dirofilaria_immitis_chr'))
data_m <- data_m %>% filter(str_detect(V1, 'dirofilaria_immitis_chr'))


plot_m <- ggplot(data_m, aes(V2, V5/median(V5)*2, col=as.factor(V1))) + geom_point(size=0.5) + scale_colour_viridis_d() + ylim(0,4) + theme_bw() + facet_grid(.~V1, scales = "free_x")

plot_f <- ggplot(data_f, aes(V2, V5/median(V5)*2, col=as.factor(V1))) + geom_point(size=0.5) + scale_colour_viridis_d() + ylim(0,4) + theme_bw() + facet_grid(.~V1, scales = "free_x")

plot_f + plot_m + plot_layout(ncol=1, guides = "collect")

ggplot(data, aes(1:nrow(data), V5/median(V5), col=as.factor(V1))) + geom_point(size=0.5) + scale_colour_viridis_d() + ylim(0,2)

# close look at the breakpoint in coverage on the X chromosome
ggplot(data, aes(V2, V5/median(V5), col=as.factor(V1))) + geom_line(size=0.5) + scale_colour_viridis_d() + ylim(0,2) + facet_grid(V1~.) + xlim(12.775e6, 12.825e6)


ggplot(data, aes(1:nrow(data), V5/median(V5), col=as.factor(V1))) + geom_point(size=0.5) + scale_colour_viridis_d() + ylim(0,2)
```



## Mapping nanopore to genome to check X chromosome coverage
- coverage breaks from low to high between 12790000 and 12800000 - likely boundary of the PAR

```bash
bsub.py 10 nanopore_to_di "minimap2 -a -x map-ont dimmitis_WSI_1.0.fa ../GENOME_IMPROVEMENT/NANOPORE_DATA/SRR14299255_1.fastq.gz  \| samtools view -h -b \| samtools sort -o nanopore_to_di.bam"

samtools index nanopore_to_di.bam


```




## Telomeres

fastaq search_for_seq dimmitis_WSI_1.0.fa telomere_dimer.pos ttaggcttaggc

cat telomere_dimer.pos | awk '{print $1,$2,$2+1}' OFS="\t" > telomere_dimer.bed

bedtools makewindows -g ${REF_DIR}/REF.genome -w 10000 > genome.10000bp.bed

bedtools coverage -b telomere_dimer.bed -a REF.genome.10000bp.bed  | awk '{print $1, $2, $3, $5}' OFS="\t" | bedtools sort > telomere_coverage.bedgraph
