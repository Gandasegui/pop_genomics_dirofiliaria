# Generate an ALL SITES variant set for running pixy properly

### Some enviromental variables and modules

```bash
#working dir
WORKING_DIR=/lustre/scratch118/infgen/team333/jg34/POP_Diro

# load gatk
module load gatk/4.1.4.1

#vcftools
module load vcftools/0.1.16-c4

# also need htslib for tabix
module load common-apps/htslib/1.9.229

# create bam list using full path to bams - this allows bams to be anywhere
ls ${WORKING_DIR}/03_MAP/*.merged.bam > ${WORKING_DIR}/04_VARIANTS/bam.list
BAM_LIST=${WORKING_DIR}/04_VARIANTS/bam.list

#The ref
REFERENCE=${WORKING_DIR}/01_REF/dimmitis_WSI_2.2.fa
REF_DIR=${WORKING_DIR}/01_REF

# make a sequences list to allow splitting jobs per scaffold/contig
grep ">" ${REFERENCE} | sed -e 's/>//g' > ${WORKING_DIR}/04_VARIANTS/sequences.list
SEQUENCE=${WORKING_DIR}/04_VARIANTS/sequences.list

# make a new directory for the merged GVCFS
mkdir ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED_ALLSITES
cd ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED_ALLSITES

# make a list of GVCFs to be merged
ls -1 ${WORKING_DIR}/04_VARIANTS/GVCFS/*complete/*gz > ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED_ALLSITES/gvcf.list
GVCF_LIST=${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED_ALLSITES/gvcf.list
```

#### Setup the run files

```bash
n=1
while read SEQUENCE; do
     echo -e "gatk CombineGVCFs -R ${REFERENCE} --intervals ${SEQUENCE} \\" > ${n}.run_merge_gvcfs.tmp.${SEQUENCE}
     while read SAMPLE; do
          echo -e "--variant ${SAMPLE} \\" >> ${n}.run_merge_gvcfs.tmp.${SEQUENCE};
     done < ${GVCF_LIST}
     echo -e "--output ${SEQUENCE}.cohort.g.vcf.gz" >> ${n}.run_merge_gvcfs.tmp.${SEQUENCE};
     let "n+=1";
done < ${WORKING_DIR}/04_VARIANTS/sequences.list
chmod a+x *.run_merge_gvcfs.tmp.*

# run
for i in *.run_merge_gvcfs.tmp.*; do
     bsub.py --queue long --threads 4 10 merge_vcfs "./${i}";
done
```

#### Split each chromosome up into separate jobs, and run genotyping on each individually.

```bash
n=1
while read SEQUENCE; do
     echo -e "gatk GenotypeGVCFs \
     -R ${REFERENCE} \
     -V ${SEQUENCE}.cohort.g.vcf.gz \
     --intervals ${SEQUENCE} \
     --all-sites \
     --heterozygosity 0.015 \
     --indel-heterozygosity 0.01 \
     --annotation DepthPerAlleleBySample --annotation Coverage --annotation ExcessHet --annotation FisherStrand --annotation MappingQualityRankSumTest --annotation StrandOddsRatio --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation DepthPerSampleHC --annotation QualByDepth \
     -O ${n}.${SEQUENCE}.cohort.vcf.gz" > run_hc_genotype.${SEQUENCE}.tmp.job_${n};
     let "n+=1";
done < ${WORKING_DIR}/04_VARIANTS/sequences.list

chmod a+x run_hc_genotype*

mkdir LOGFILES

# setup job conditions
JOBS=$( ls -1 run_hc_* | wc -l )
ID="U$(date +%s)"

# run
bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 4 -M20000 -J GATK_HC_GENOTYPE_${ID}_[1-$JOBS] -e LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].e -o LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].o "./run_hc_*\$LSB_JOBINDEX"
```

#### Bring the files together

```bash
# make list of vcfs
ls -1 *.cohort.vcf.gz | sort -n > vcf_files.list

# merge them
vcf-concat --files vcf_files.list > dirofilaria_immitis.cohort.allsites.vcf;
     bgzip dirofilaria_immitis.cohort.allsites.vcf;
     tabix -p vcf dirofilaria_immitis.cohort.allsites.vcf.gz

# clean up
rm run*
rm ^[0-9]*
rm *.g.vcf.gz*
```

#### Now, let's filter nuclear variants and invariants

```bash
VCF=${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED_ALLSITES/dirofilaria_immitis.cohort.allsites.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb

vcftools --gzvcf dirofilaria_immitis.cohort.allsites.vcf.gz --remove-indels
#After filtering, kept 32 out of 32 Individuals
#After filtering, kept 85744962 out of a possible 86372104 Sites

#select nuclear invariants
bsub.py 1 select_nuclearINVARIANTs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include NO_VARIATION \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf"

vcftools --vcf dirofilaria_immitis.cohort.allsites.nuclearINVARIANTs.vcf --remove-indels
#After filtering, kept 32 out of 32 Individuals
#After filtering, kept 85270107 out of a possible 85310047 Sites
```

Filtering has an strong effect int the invarins and does not work so well

Some guides reocommned not to apply pop genomics filters for invariants

#### In any way, lets merge the nuclear invariants with nuclear SNPs previously filtered

```bash
#Keep only the 19 samples of the nuclear dataset and then merge with the SNPs-INDELs
bsub.py 1 filter_nuclear_INVARITANTvcftools "vcftools --vcf ${VCF%.vcf.gz}.nuclearINVARIANTs.vcf \
--keep ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED/FILTER/nuclear_samplelist.keep \
--recode --out nuclearINVARIANTS_19samples"
#After filtering, kept 19 out of 32 Individuals
#After filtering, kept 85310047 out of a possible 85310047 Sites

#And merge with the already filtered variants
bsub.py --done "filter_nuclear_INVARITANTvcftools" 1 merge_nuclear_VARIANTsandINVARIANTs "gatk MergeVcfs \
--INPUT ${WORKING_DIR}/04_VARIANTS/FINAL_SETS/nuclear_samples3x_missing0.8.recode.vcf \
--INPUT nuclearINVARIANTS_19samples.recode.vcf \
--OUTPUT nuclearVARIANTsandINVARIANTs_19samples.recode.vcf"

#Also, we will slectct only the chrX to chr4, avoiding the scaffolds
bsub.py 1 filter_CHR "vcftools \
--vcf nuclearVARIANTsandINVARIANTs_19samples.recode.vcf \
--chr dirofilaria_immitis_chrX \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--remove-indels \
--recode --out nuclearSNPssandINVARIANTs.chrxto4"
#After filtering, kept 19 out of 19 Individuals
#After filtering, kept 85335439 out of a possible 85732982 Sites

bgzip nuclearSNPssandINVARIANTs.chrxto4.recode.vcf;
tabix -p vcf nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz
mv nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.* ../FINAL_SETS/
```

## Now, we run pixy using that file

```bash
#Create and set the environment
conda create --name pixy
conda init --all
conda activate pixy

#Installing pixy
conda install -c conda-forge pixy
conda install -c bioconda htslib

#Few environemntal variables
WORKING_DIR=/lustre/scratch118/infgen/team333/jg34/POP_Diro
VCF=${WORKING_DIR}/04_VARIANTS/FINAL_SETS/nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz

mkdir ${WORKING_DIR}/05_ANALYSIS
mkdir ${WORKING_DIR}/05_ANALYSIS/PIXY
cd ${WORKING_DIR}/05_ANALYSIS/PIXY

#Generate two population files
#country_pop.list
'AUS_NSW_AD_001	AUS
AUS_NSW_AD_002	AUS
AUS_NSW_AD_003	AUS
AUS_NSW_AD_004	AUS
AUS_NSW_AD_005	AUS
ITL_PAV_ADF_001	ITL
USA_ARK_MFP_001	USA
USA_GEO_ADF_001	USA
USA_GEO_MFP_001	USA
USA_GEO_MFP_002	USA
USA_GEO_MFS_001	USA
USA_GEO_MFS_002	USA
USA_ILL_MFS_001	USA
USA_LOU_MFP_001	USA
USA_LOU_MFS_001	USA
USA_MCH_MFP_001	USA
USA_MIP_MFS_001	USA
USA_TEN_MFP_001	USA
USA_TEX_MFP_001	USA'
#sampletype_pop.list
'AUS_NSW_AD_001	single
AUS_NSW_AD_002	single
AUS_NSW_AD_003	single
AUS_NSW_AD_004	single
AUS_NSW_AD_005	single
ITL_PAV_ADF_001	single
USA_ARK_MFP_001	pooled
USA_GEO_ADF_001	single
USA_GEO_MFP_001	pooled
USA_GEO_MFP_002	pooled
USA_GEO_MFS_001	single
USA_GEO_MFS_002	single
USA_ILL_MFS_001	single
USA_LOU_MFP_001	pooled
USA_LOU_MFS_001	single
USA_MCH_MFP_001	pooled
USA_MIP_MFS_001	single
USA_TEN_MFP_001	pooled
USA_TEX_MFP_001	pooled'

#Submit the two jobs to get the diversity per country and sampletype
bsub.py --queue long --threads 20 20 pixy_country \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations country_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix country"

bsub.py --queue long --threads 20 20 pixy_sampletype \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations sampletype_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix sampletype"
```

#### Download that files into R an proceed to estimate the values and generate plots

```R
library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggridges)
```

### Nucleotide diversity (pi)

```R
##### PI #####

# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("data/country_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS")
pi_data_USA <- pi_data %>%
  filter(pop=="USA")
pi_data_ITL <- pi_data %>%
  filter(pop=="ITL")

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

# get position for vertical lines used in plot to delineate the linkage groups
pi_data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

'
# A tibble: 5 × 2
chromosome   max
<chr>      <int>
  1 chr1         440
2 chr2         592
3 chr3         744
4 chr4         885
5 chrX         283
'

#LEt's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 2 × 2
chr_type   median
<chr>       <dbl>
1 autosome 0.000269
2 sexchr   0.000133
'
# 0.000133 / 0.000269 = 0.4944238 (far off 0.75 expected of diversity on sex chromosome relative to autosome)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 6 × 3
# Groups:   pop [3]
pop   chr_type     median
<chr> <chr>         <dbl>
1 AUS   autosome 0.0000391 
2 AUS   sexchr   0.00000376
3 ITL   autosome 0.0000864 
4 ITL   sexchr   0.0000322 
5 USA   autosome 0.000736  
6 USA   sexchr   0.000426  
'
# plot 1 - genome wide plots per population
plot_1 <- ggplot(pi_data, aes(position*100000, avg_pi, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(pop~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Nucleotide Diversity (Pi)")


# plot 2 - density plots of pi per group
plot_2 <- ggplot(pi_data, aes(avg_pi, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(pop~.) +
  xlim(0, 0.005) +
  scale_x_continuous(breaks = c(0, 0.0025, 0.005)) +
  scale_fill_npg() +
  labs(x="Nucleotide Diversity (Pi)", y="Density")

# bring it together
plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("FigS7_plots_genomewide_and_density_Pi.png", width=9, height=6)

#Now a boxplot of the pi value per population

scale_colour_javi <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('dodgerblue4', 'turquoise4', 'darkgoldenrod1'), 
      c('AUS', 'ITL', 'USA')), 
    ...
  )
}

boxplot_pi <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_javi ()+
  ylim(0, 0.003)

boxplot_pi

ggsave("Fig2c_plots_boxplot_pop_Pi.png", width=4, height=4)

#let's check the statistical significance of the differences between pi values
#Let's explore the normality
pi_data %>% filter(chromosome != 'chrX') %>%
ggplot(., aes(x = avg_pi, colour = pop)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ pop)
#Confirm the absence of normality by shapiro
pi_data_AUS <- pi_data_AUS %>% filter(chromosome != 'chrX')
AUS_shapiro <- shapiro.test(pi_data_AUS$avg_pi)
print(AUS_shapiro)#W = 0.69062, p-value < 2.2e-16
pi_data_ITL <- pi_data_ITL %>% filter(chromosome != 'chrX')
ITL_shapiro <- shapiro.test(pi_data_ITL$avg_pi)
print(ITL_shapiro)#W = 0.69994, p-value < 2.2e-16
pi_data_USA <- pi_data_USA %>% filter(chromosome != 'chrX')
USA_shapiro <- shapiro.test(pi_data_USA$avg_pi)
print(USA_shapiro)#W = 0.90598, p-value < 2.2e-16

#Willcoxson test for everyone
wilcox.test(pi_data_AUS$avg_pi, pi_data_ITL$avg_pi)
#W = 152048, p-value = 1.16e-06
wilcox.test(pi_data_AUS$avg_pi, pi_data_USA$avg_pi)
#W = 41288, p-value < 2.2e-16
wilcox.test(pi_data_ITL$avg_pi, pi_data_USA$avg_pi)
#W = 50720, p-value < 2.2e-16

#let's generate some dataframe for the estatistics
# subset
pi_data_AUS <- pi_data %>%
  filter(pop=="AUS") %>%
  filter(chromosome != 'chrX')

pi_data_USA <- pi_data %>%
  filter(pop=="USA") %>%
  filter(chromosome != 'chrX') 

pi_data_ITL <- pi_data %>%
  filter(pop=="ITL") %>%
  filter(chromosome != 'chrX')


```

# Dxy and Fst

```R
# load data
dxy_data <- read.table("data/country_dxy.txt", header=T)
fst_data <- read.table("data/country_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)
data$chromosome <- str_remove(data$chromosome, 'dirofilaria_immitis_')

# add numbering to help with plotting
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())

# add sex chromosome information
data <- data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))


# summarise median data for Fst and Dxy
data %>%
  group_by(comparison,data_type) %>%
  summarise(median = median(value, na.rm = TRUE))

'
# A tibble: 6 × 3
# Groups:   comparison [3]
comparison data_type   median
<chr>      <chr>        <dbl>
1 AUS_v_ITL  Dxy       0.000608
2 AUS_v_ITL  Fst       0.876   
3 AUS_v_USA  Dxy       0.000643
4 AUS_v_USA  Fst       0.327   
5 ITL_v_USA  Dxy       0.000666
6 ITL_v_USA  Fst       0.161
'

#Plotting Dxy

# plot 1 - genome wide plots per comparison
plot_1 <- data %>%
  filter(data_type =='Dxy')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Dxy")


# plot 2 - density plots of dxy per group
plot_2 <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Dxy", y="Density")

# bring it together
dxy_plot <- plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_dxy.png", width=9, height=6)

#some additional plots
boxplot_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(comparison, value, col=comparison)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="black") +
  labs(x = "Population" , y = "Dxy", colour = "Population") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

density_dxy <- data %>%
  filter(data_type =='Dxy') %>% 
  ggplot(., aes(x=value, group = comparison, fill=comparison)) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Dxy" , y = "density", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

scatter_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ITL, y = AUS_v_USA, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ITL" , y = "AUS_v_USA", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
scatter_dxy_b <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  ggplot(., aes(x=AUS_v_ITL, y = ITL_v_USA, col=chromosome)) +
  geom_point(alpha=.4) +
  labs(x = "AUS_v_ITL" , y = "ITL_v_USA", colour = "Comparison") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

genome_pos_dxy_a <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL - AUS_v_USA,
         "yy" = AUS_v_ITL - ITL_v_USA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL - AUS vs USA")

genome_pos_dxy_b <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL - AUS_v_USA,
         "yy" = AUS_v_ITL - ITL_v_USA) %>%
  ggplot(., aes(position*100000, yy, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL - ITL vs USA")

genome_pos_dxy_aa <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL / AUS_v_USA,
         "yy" = AUS_v_ITL / ITL_v_USA) %>%
  ggplot(., aes(position*100000, xx, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL / AUS vs USA")

genome_pos_dxy_bb <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL / AUS_v_USA,
         "yy" = AUS_v_ITL / ITL_v_USA) %>%
  ggplot(., aes(position*100000, yy, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS vs ITL / ITL vs USA")

genome_pos_dxy_cc <-  data %>%
  filter(data_type == 'Dxy') %>%
  as.data.frame() %>%
  select('position', 'comparison', 'value', 'chromosome')%>%
  as_tibble() %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate("xx" = AUS_v_ITL / AUS_v_USA,
         "yy" = AUS_v_ITL / ITL_v_USA,
         "zz" = AUS_v_USA / ITL_v_USA) %>%
  ggplot(., aes(position*100000, zz, col=chromosome)) +
  geom_point(size=1) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Position", y="AUS_v_USA / ITL vs USA")
  
ggarrange(scatter_dxy_a, scatter_dxy_b, common.legend = T)
ggarrange(genome_pos_dxy_a, genome_pos_dxy_b, common.legend = T, ncol = 1)
ggarrange(genome_pos_dxy_aa, genome_pos_dxy_bb, genome_pos_dxy_cc,
          common.legend = T, ncol = 1)

ggsave("nosequecosa.png", width=9, height=6)
#Trying to plot lines
#getting chr positions
chr <- c('X', '1', '2', '3', '4')

for (i in chr) {
  x <- data %>%
    filter(chromosome == paste0('chr', i))
  print(quantile(x$position))
  }

data$comparison <- str_replace_all(data$comparison, '_', ' ')
data$comparison <- str_replace_all(data$comparison, 'v', 'vs')

dxy_lineplot <- data %>%
  filter(data_type=='Dxy')%>%
  ggplot(., aes(position*100000, value, col=comparison)) +
  geom_point(size=0.2, alpha = 0.1) +
  geom_vline(xintercept=c(281*100000,
                          438*100000,
                          590*100000,
                          740*100000),
                          size=1, linetype="dashed", col='grey41')+
  annotate(geom="text", x=140*100000, y=0.004, label="ChrX")+
  annotate(geom="text", x=359.5*100000, y=0.004, label="Chr1")+
  annotate(geom="text", x=545*100000, y=0.004, label="Chr2")+
  annotate(geom="text", x=665*100000, y=0.004, label="Chr3")+
  annotate(geom="text", x=810.5*100000, y=0.004, label="Chr4")+
  geom_line() +
  scale_color_brewer(type = 'div', palette = 'Accent') +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Dxy")

dxy_lineplot
# Plotting Fst

# plot 1 - genome wide plots per population
plot_1 <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(., aes(position*100000, value, col=chromosome)) +
  geom_point(size=1) +
  facet_grid(comparison~.) +
  scale_color_tron() +
  theme_bw() +
  theme(legend.position="top", 
        legend.title = element_blank()) +
  labs(x="Genomic Position", y="Fst")


# plot 2 - density plots of Fst per group
plot_2 <- data %>%
  filter(data_type=='Fst')%>%
  ggplot(aes(value, chr_type, fill=chr_type), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +
  theme_bw() + theme(legend.position = "none", axis.text.y=element_blank()) +
  facet_grid(comparison~.) +
  scale_fill_npg() +
  labs(x="Fst", y="Density")

# bring it together
plot_1 + plot_2 +  plot_layout(widths = c(5, 1))
ggsave("plots_genomewide_and_density_fst.png", width=9, height=6)
```

### Now let's ecplore the american populations


```bash
#Let's select only the USA samples
cd ${WORKING_DIR}/04_VARIANTS/FINAL_SETS/

vcftools --gzvcf nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz \
--keep USA_samplelist.keep \
--recode --out USA_nuclearSNPssandINVARIANTs.chrxto4

bgzip USA_nuclearSNPssandINVARIANTs.chrxto4.recode.vcf;
tabix -p vcf USA_nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz

#Now, let's run pixy
VCF=${WORKING_DIR}/04_VARIANTS/FINAL_SETS/USA_nuclearSNPssandINVARIANTs.chrxto4.recode.vcf.gz

#usa_pop.list
USA_ARK_MFP_001	ARK
USA_GEO_ADF_001	GEO
USA_GEO_MFP_001	GEO
USA_GEO_MFP_002	GEO
USA_GEO_MFS_001	GEO
USA_GEO_MFS_002	GEO
USA_ILL_MFS_001	ILL
USA_LOU_MFP_001	LOU
USA_LOU_MFS_001	LOU
USA_MCH_MFP_001	MCH
USA_MIP_MFS_001	MIP
USA_TEN_MFP_001	TEN
USA_TEX_MFP_001	TEX

#usa_sampletype.list
USA_ARK_MFP_001	pooled
USA_GEO_ADF_001	single
USA_GEO_MFP_001	pooled
USA_GEO_MFP_002	pooled
USA_GEO_MFS_001	single
USA_GEO_MFS_002	single
USA_ILL_MFS_001	single
USA_LOU_MFP_001	pooled
USA_LOU_MFS_001	single
USA_MCH_MFP_001	pooled
USA_MIP_MFS_001	single
USA_TEN_MFP_001	pooled
USA_TEX_MFP_001	pooled

#and submitting jobs
bsub.py --queue long --threads 20 20 usa_pop \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations usa_pop.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix usa_pop"

bsub.py --queue long --threads 20 20 usa_sampletype \
"pixy --stats pi fst dxy \
--vcf ${VCF} \
--populations usa_sampletype.list \
--window_size 100000 \
--n_cores 20 \
--output_prefix usa_sampletype"
```

### Now, let's import into R and generate some plots

```R
##### PI #####

# Per population
# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("data/usa_pop_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#LEt's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))

# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 2 × 2
chr_type   median
<chr>       <dbl>
1 autosome 0.000697
2 sexchr   0.000368
'
# 0.000133 / 0.000269 = 0.527977 (far off 0.75 expected of diversity on sex chromosome relative to autosome)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 16 × 3
# Groups:   pop [8]
pop   chr_type   median
<chr> <chr>       <dbl>
  1 ARK   autosome 0.00134 
2 ARK   sexchr   0.000651
3 GEO   autosome 0.000781
4 GEO   sexchr   0.000457
5 ILL   autosome 0.000261
6 ILL   sexchr   0.000130
7 LOU   autosome 0.000869
8 LOU   sexchr   0.000451
9 MCH   autosome 0.00117 
10 MCH   sexchr   0.000581
11 MIP   autosome 0.000235
12 MIP   sexchr   0.000146
13 TEN   autosome 0.00126 
14 TEN   sexchr   0.000610
15 TEX   autosome 0.000681
16 TEX   sexchr   0.000277
'

#Now a boxplot of the pi value per population
pi_data$pop <- factor(pi_data$pop, 
                           levels = c('MCH', 'ILL', 'TEN', 'ARK', 
                                      'GEO', 'MIP', 'TEX', 'LOU'))

scale_colour_javier_PCA <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('royalblue1', 'blue3',
        'turquoise3',
        'green',
        'lightgoldenrod1', 'darkgoldenrod1', 'orange', 'orange2',
        'orange3', 'tomato1', 'red2', 'darkred'), 
      c('QUE', 'NSW', 
        'PAV',
        'SIC',
        'MCH', 'ILL', 'TEN', 'ARK', 
        'GEO', 'MIP', 'TEX', 'LOU')), 
    ...
  )
}

boxplot_pi_usa_pop <- ggplot(pi_data, aes(pop, avg_pi, col=pop)) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_boxplot(fill=NA, col="grey15", linewidth = 1) +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_npg() +
  ylim(0, 0.003) + scale_colour_javier_PCA()

#Per sample type
# get nucleotide diversity (pi) data from pixy output
pi_data <- read.table("data/usa_sampletype_pi.txt", header=T)
pi_data$chromosome <- str_remove(pi_data$chromosome, 'dirofilaria_immitis_')

# filter pi data to remove small scaffolds not in the linkage groups, and to number the rows per group to help with plotting
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

#LEt's add the chr type variable
pi_data <- pi_data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "X"), "sexchr", "autosome"))
# calculate the median Pi and checking the ratio of sex-to-autosome diversity. Should be about 0.75, as Trichuris is XX/XY
pi_data_sex_median <- pi_data %>%
  group_by(chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 2 × 2
chr_type   median
<chr>       <dbl>
1 autosome 0.000704
2 sexchr   0.000378
'
# 0.000378 / 0.000704 = 0.5369318 (far off 0.75 expected of diversity on sex chromosome relative to autosome)

pi_data_pop_sex_median <- pi_data %>%
  group_by(pop, chr_type) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))
'
# A tibble: 4 × 3
# Groups:   pop [2]
  pop    chr_type   median
  <chr>  <chr>       <dbl>
1 pooled autosome 0.000881
2 pooled sexchr   0.000480
3 single autosome 0.000539
4 single sexchr   0.000315
'
#Now a boxplot of the pi value per population
boxplot_pi_usa_sampletype <- ggplot(pi_data, aes(pop, avg_pi, col=pop, shape = pop)) +
  geom_jitter(size = 2, alpha = 0.5) +
  geom_boxplot(fill=NA, col="grey15", linewidth = 1) +
  labs(x = "Population" , y = "Nucleotide diversity (Pi)") +
  theme_bw() +
  scale_color_tron()+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0, 0.003)

#and arranging
ggarrange(boxplot_pi_usa_pop, boxplot_pi_usa_sampletype, labels = c('a', 'b'),
          ncol = 2, widths =  c(1.5, 1))
ggsave("FigS8_plots_diversity_USA.png", width=9, height=5.5)
```
