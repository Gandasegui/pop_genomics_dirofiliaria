# R script to carry out the PCA

```R
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(ggsci)
library(ggpubr)
library(reshape2)
library(viridis)
library(vcfR)
library(factoextra)

#Preparing the data for plotin
scale_colour_javier_mitoPCA <- function(...){
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

#PCA on nuclear variants using genotypes
snpgdsClose(genofile)
vcf.in <- "data/nuclear_samples3x_missing0.8.chr1to4.recode.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "nucDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)


pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"
metadata <- samples %>% separate(name,c("country","population","sampletype","sampleID"))
metadata$sampletype <- gsub('ADF', 'AD', metadata$sampletype)

data <- data.frame(sample.id = pca$sample.id,
                   EV1 = pca$eigenvect[,1],  
                   EV2 = pca$eigenvect[,2],
                   EV3 = pca$eigenvect[,3],
                   EV4 = pca$eigenvect[,4],
                   EV5 = pca$eigenvect[,5],
                   EV6 = pca$eigenvect[,6],
                   COUNTRY = metadata$country,
                   POPULATION = metadata$population,
                   SAMPLETYPE = metadata$sampletype,
                   SAMPLEID = metadata$sampleID,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c('QUE', 'NSW', 
                                     'PAV',
                                     'SIC',
                                     'MCH', 'ILL', 'TEN', 'ARK', 
                                     'GEO', 'MIP', 'TEX', 'LOU'))

# Plot PC1 vs PC2
nuc_PC1 <- ggplot(data,aes(EV1, EV2, col = POPULATION, label = POPULATION)) +
  geom_point(alpha = 0.8, size = 3) +
  theme_bw() +
  scale_colour_javier_mitoPCA () +
  labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
       y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"))

#Plot PC3 vs PC4
nuc_PC3 <- ggplot(data,aes(EV3, EV4, col = POPULATION, label = POPULATION)) +
  geom_point(alpha = 0.8, size = 3) + 
  theme_bw() +
  scale_colour_javier_mitoPCA () +
  labs(x = paste0("PC3 variance: ",round(pca$varprop[3]*100,digits=2),"%"),
       y = paste0("PC4 variance: ",round(pca$varprop[4]*100,digits=2),"%"))
```
We have obeserved some substructure that can be due to seq noise

We will findout what is going on

### Let's use SNPs freqfr PCA and extract tecnical information:

-Sample type (pooled or single)

-Depth

-Missingness

-Heterozygosity (inbreeding coeficient)

```R
#Now we load the vcf file, just to understand the structure
#We use the file with no indels
#Read two times the vcf file, first for the columns names, second for the data
nuc_vcf<-readLines("data/nuclear_samples3x_missing0.8.chr1to4.NOindels.recode.vcf")
nuc_vcf_data<-read.table("data/nuclear_samples3x_missing0.8.chr1to4.NOindels.recode.vcf", 
                         stringsAsFactors = FALSE)
# filter for the columns names
nuc_vcf<-nuc_vcf[-(grep("CHROM",nuc_vcf)+1):-(length(nuc_vcf))]
vcf_names<-unlist(strsplit(nuc_vcf[length(nuc_vcf)],"\t"))
names(nuc_vcf_data)<-vcf_names

#To obtain parameters, let's use the vcfR package
#Let's extrct the depth for each sample
vcf_file <- system.file("extdata", "data/nuclear_samples3x_missing0.8.chr1to4.NOindels.recode.vcf", package = "pinfsc50")
vcf <- read.vcfR("data/nuclear_samples3x_missing0.8.chr1to4.NOindels.recode.vcf", verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
dp <- as.data.frame(dp)
dp[] <- lapply(dp, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

#Now, we have created a df with the depth for each sample and each SNP
#Let's estimate de mean depth and create a dataframe to merge later on
file_names <- colnames(dp)

#The , I create an empty df
dp_stats <- data_frame(name = character(), dp_mean = double(), dp_sd = double())

for (i in 1:19) {
  #first, I create the name
  name <- file_names[i]
  x <- dp[, i]
  # reading the table
  m <- mean(x, na.rm = T)
  sd <- sd(x, na.rm = T)
  #now, I add the row to the df
  dp_stats <- dp_stats %>% add_row('name' = name, 'dp_mean' = m, 'dp_sd' = sd)
  #print(nuc_cov)
  print(dp_stats)
}

#Let's check if missingess is affecting PCA
#Although we estimated the missingness before, we should estimate it again
#with the new set of SNPs -- (see bash scritp)

#vcftools --vcf nuclear_samples3x_missing0.8.chr1to4.recode.vcf \
#--out nuclear_final --missing-indv

#First we read the table with missingness info
missing <- read.table('data/nuclear_final.imiss', header = T) %>%
  select(INDV, F_MISS)
colnames(missing) <- c('name', 'miss')

#Now, we estimate the herterozygosity

#vcftools --vcf nuclear_samples3x_missing0.8.chr1to4.recode.vcf \
#--het --out nuclear_final

#And read the table
het <- read.table('data/nuclear_final.het', header = T) %>%
  select(INDV, 'F')
colnames(het) <- c('name', 'het')

#Let's run the PCA using SNP frequencies
#First, we will extract the AD, which refers to 
ad <- extract.gt(vcf, element='AD', as.numeric=F)
ad <- as.data.frame(ad)
file_names <- colnames(ad)
#Let's do it in a loop
#First we create an empty df with the same rownames
allele_freq <- data.frame(matrix(, nrow=216455, ncol=0))
rownames(allele_freq) <- rownames(ad)

#And now I will do it in a loop
for (i in 1:19) {
  #first, I create the name
  name <- file_names[i]
  #The, I play with the data to obtain a new freq column
  x <- as_vector(ad[, i])
  z<- as_tibble(str_split_fixed(x, ",", 2))
  z[] <- lapply(z, function(x) {
    if(is.character(x)) as.numeric(as.character(x)) else x
  })
  z<- z %>% mutate(V3 = V2/(V1+V2)*100)
  #Add col names based on 'name'
  colnames(z) <- c(paste0(name, '_ref'), paste0(name, '_alt'), paste0(name, '_freq'))
  #now, I add the row to the df
  allele_freq <-allele_freq %>% cbind(z)
}

#Now we select those frequencies to perfomr pca
freqs <- allele_freq %>%
  select(., contains('_freq'))

freq.pca <- allele_freq %>%
  select(., contains('_freq')) %>%
  na.omit() %>% #after omiting NAs, 126728 variants are retained
  prcomp(., scale = TRUE)

#Let's see the screeplot
fviz_eig(freq.pca)
summary(freq.pca)

freq.pca.eigenvect <- as.data.frame(freq.pca$rotation)
rownames(freq.pca.eigenvect) <- gsub('_freq', '', rownames(freq.pca.eigenvect))
freq.pca.eigenvect <- rownames_to_column(freq.pca.eigenvect, var = 'name')
samples <- as.data.frame(freq.pca.eigenvect$name)
colnames(samples) <- "name"
samples <- left_join(samples, missing, by = 'name')
samples <- left_join(samples, dp_stats, by = 'name')
samples <- left_join(samples, het, by = 'name')
metadata <- samples %>% separate(name,c("country","population","sampletype","sampleID"))
freq.pca.eigenvect <- cbind(freq.pca.eigenvect, metadata)
freq.pca.eigenvect <- freq.pca.eigenvect %>% mutate(sampletype_2 = recode(sampletype,
                                                      'AD' = 'Single',
                                                      'ADF' = 'Single',
                                                      'MFP' = 'Pooled',
                                                      'MFS' = 'Single'))


#Generating levels for population and super-population variables
freq.pca.eigenvect$population <- factor(freq.pca.eigenvect$population, 
                          levels = c('QUE', 'NSW', 
                                     'PAV',
                                     'SIC',
                                     'MCH', 'ILL', 'TEN', 'ARK', 
                                     'GEO', 'MIP', 'TEX', 'LOU'))

#Adding new variables
colnames(freq.pca.eigenvect) <- str_replace(colnames(freq.pca.eigenvect), 'population', 'Population')
colnames(freq.pca.eigenvect) <- str_replace(colnames(freq.pca.eigenvect), 'sampletype_2', 'Sample_type')
colnames(freq.pca.eigenvect) <- str_replace(colnames(freq.pca.eigenvect), 'dp_mean', 'Depth')
colnames(freq.pca.eigenvect) <- str_replace(colnames(freq.pca.eigenvect), 'miss', 'Missingness')
colnames(freq.pca.eigenvect) <- str_replace(colnames(freq.pca.eigenvect), 'het', 'Heterozygosity')

#Let's plot again
allele_freq_PCA <- ggplot(freq.pca.eigenvect, aes(PC1, PC2, col = Population, label = name)) +
  geom_point(size=4) +
  theme_bw() + scale_colour_javier_mitoPCA() +
  theme(legend.position="none")+
  labs(x = paste0("PC1 variance: 54.38%"),
       y = paste0("PC2 variance: 15.03%"))
       
allele_freq_PCA
ggsave('Figures/Fig2b_PCA_nuc.png')

#First we will plot the the zoom on USA samples
#With sampletype 
stype <- ggplot(freq.pca.eigenvect,aes(PC1, PC2, col = Population, shape = Sample_type)) +
  geom_point(size=4) +
  theme_bw() + 
  scale_colour_javier_mitoPCA () +
  labs(x = "PC1",
       y = "PC2") +
  xlim(0.175, 0.30) +
  ylim(0.1, 0.2)
#With sample type + DP
dp <- ggplot(freq.pca.eigenvect,aes(PC1, PC2, col = log(Depth), shape = Sample_type)) +
  geom_point(size=4) +
  theme_bw() + 
  scale_colour_viridis() +
  labs(x = "PC1",
       y = "PC2") +
  xlim(0.175, 0.30) +
  ylim(0.1, 0.2)
#With sample type + missingess
miss <- ggplot(freq.pca.eigenvect,aes(PC1, PC2, col = Missingness, shape = Sample_type)) +
  geom_point(size=4) +
  theme_bw() + 
  scale_colour_viridis() +
  labs(x = "PC1",
       y = "PC2") +
  xlim(0.175, 0.30) +
  ylim(0.1, 0.2)
#With sample type + heterozygosity
hetero <- ggplot(freq.pca.eigenvect,aes(PC1, PC2, col = Heterozygosity, shape = Sample_type)) +
  geom_point(size=4) +
  theme_bw() + 
  scale_colour_viridis() +
  labs(x = "PC1",
       y = "PC2") +
  xlim(0.175, 0.30) +
  ylim(0.1, 0.2)

ggarrange(stype, dp, miss, hetero, nrow = 2, common.legend = F, 
          labels=c("a",'b', 'c', 'd'), vjust = 1, ncol = 2)
ggsave("Figures/FigS6_PCA_DP_MISS_HET.jpg", height = 8, width = 9) 
```

### Now let's run the PCA using genotypes in the mito and Wb variants

```R
#snpgdsClose(genofile)
vcf.in <- "data/mito_samples3x_missing0.8.recode.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "mtDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)

pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"
metadata <- samples %>% separate(name,c("country","population","sampletype","sampleID"))

data <- data.frame(sample.id = pca$sample.id,
                   EV1 = pca$eigenvect[,1],  
                   EV2 = pca$eigenvect[,2],
                   EV3 = pca$eigenvect[,3],
                   EV4 = pca$eigenvect[,4],
                   EV5 = pca$eigenvect[,5],
                   EV6 = pca$eigenvect[,6],
                   COUNTRY = metadata$country,
                   POPULATION = metadata$population,
                   SAMPLETYPE = metadata$sampletype,
                   SAMPLEID = metadata$sampleID,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c('QUE', 'NSW', 
                                     'PAV',
                                     'SIC',
                                     'MCH', 'ILL', 'TEN', 'ARK', 
                                     'GEO', 'MIP', 'TEX', 'LOU'))

# Plot PC1 vs PC2
mito_PC1 <- ggplot(data,aes(EV1, EV2, col = POPULATION, label = POPULATION)) +
  geom_point(alpha = 0.8, size = 3) + 
  theme_bw() +
  scale_colour_javier_mitoPCA () +
  labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
       y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"))

##############################

#PCA on wb variants
snpgdsClose(genofile)
vcf.in <- "data/wb_samples3x_missing0.7.recode.vcf"
gds<-snpgdsVCF2GDS(vcf.in, "wbDNA.gds", method="biallelic.only")
genofile <- snpgdsOpen(gds)

pca <-snpgdsPCA(genofile, num.thread=2, autosome.only = F)
samples <- as.data.frame(pca$sample.id)
colnames(samples) <- "name"
metadata <- samples %>% separate(name,c("country","population","sampletype","sampleID"))

data <- data.frame(sample.id = pca$sample.id,
                   EV1 = pca$eigenvect[,1],  
                   EV2 = pca$eigenvect[,2],
                   EV3 = pca$eigenvect[,3],
                   EV4 = pca$eigenvect[,4],
                   EV5 = pca$eigenvect[,5],
                   EV6 = pca$eigenvect[,6],
                   COUNTRY = metadata$country,
                   POPULATION = metadata$population,
                   SAMPLETYPE = metadata$sampletype,
                   SAMPLEID = metadata$sampleID,
                   stringsAsFactors = FALSE)

#Generating levels for population and super-population variables
data$POPULATION <- factor(data$POPULATION, 
                          levels = c('QUE', 'NSW', 
                                     'PAV',
                                     'SIC',
                                     'MCH', 'ILL', 'TEN', 'ARK', 
                                     'GEO', 'MIP', 'TEX', 'LOU'))

# Plot PC1 vs PC2
wb_PC1 <- ggplot(data,aes(EV1, EV2, col = POPULATION, label = POPULATION)) +
  geom_point(alpha = 0.8, size = 3) + 
  theme_bw() +
  scale_colour_javier_mitoPCA () +
  labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
       y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"))

ggarrange(mito_PC1, wb_PC1, nrow = 1, common.legend = T, 
          labels=c("a",'b'), vjust = 1, ncol = 2)
ggsave("Figures/FigS5_PCA_mito_wb.jpg", height = 5.5, width = 9)
```
### Now let's look for variants that can allow us to differenciated the AUS samples

```R
#Read two times the vcf file, first for the columns names, second for the data
mit_vcf<-readLines("data/mito_samples3x_missing0.8.recode.vcf")
mit_vcf_data<-read.table("data/mito_samples3x_missing0.8.recode.vcf", 
                         stringsAsFactors = FALSE)

vcf <- read.vcfR("data/mito_samples3x_missing0.8.recode.vcf", verbose = FALSE)
mito_gen_AUS <- as_tibble(extract.gt(vcf)) %>% select(contains('AUS'))
# filter for the columns names
mit_vcf<-mit_vcf[-(grep("CHROM",mit_vcf)+1):-(length(mit_vcf))]
vcf_names<-unlist(strsplit(mit_vcf[length(mit_vcf)],"\t"))
names(mit_vcf_data)<-vcf_names
mit_vcf_data <- as_tibble(mit_vcf_data)
colnames(mit_vcf_data) <- str_remove(colnames(mit_vcf_data), '#')
mito_gen_AUS <-mit_vcf_data %>%
  select(POS, REF, ALT) %>%
  cbind(., mito_gen_AUS) %>%
  as_tibble()

write_csv(mito_gen_AUS, 'mito_AUS_genotypes.csv')
 ```
 The differences in the PCA are due to missing SNPs and other variants, but nothing to discriminate between pops
