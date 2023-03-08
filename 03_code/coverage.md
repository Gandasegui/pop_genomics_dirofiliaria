# Script for estimating an plotting the coverage

```bash
WORKING_DIR=/lustre/scratch118/infgen/team333/jg34/POP_Diro/
cd ${WORKING_DIR}/03_MAP/

WINDOW='100000'

for i in *.bam; do

bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"1"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.bed
bamtools header -in ${i} | grep "^@SQ" | awk -F'[:\t]' '{printf $3"\t"$5"\n"}' OFS="\t" > ${i%.bam}.chr.genome

bedtools makewindows -g ${i%.bam}.chr.genome -w ${WINDOW} > ${i%.bam}.${WINDOW}_window.bed

samtools bedcov -Q 20 ${i%.bam}.chr.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.chr.cov
samtools bedcov -Q 20 ${i%.bam}.${WINDOW}_window.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.${WINDOW}_window.cov

rm ${i%.bam}.chr.bed ${i%.bam}.${WINDOW}_window.bed ${i%.bam}.chr.genome;

done

for i in *.chr.cov; do 

printf "${i}\n" > ${i}.tmp | awk '{print $5}' OFS="\t" ${i} >> ${i}.tmp;

done

paste *.tmp > coverage_stats.summary
rm *.tmp

mkdir COV_STATS
mv *.chr.cov *_window.cov *.cov coverage_stats.summary COV_STATS/
cd COV_STATS/
```

#### Generate quantitative stats on coverage for supplementary tables etc

Extract mtDNA, Wb and nuclear (mean & stddev) data

For nuclear, we will select only the defined Chr (chrX and chr1 to chr4)

```bash
# extract mtDNA and nuclear (mean & stddev) data
for i in *.chr.cov; do
	name=${i%.chr.cov};
	nuc=$(grep -v "scaffold\|Wb\|Mt" ${i%.merged.chr.cov}.merged.100000_window.cov | datamash mean 5 sstdev 5 );
	mtDNA=$(grep "chrMtDNA" ${i} | cut -f5 );
	Wb=$(grep 'chrWb' ${i} | cut -f5 ); 
	echo -e "${name}\t${nuc}\t${mtDNA}\t${Wb}";
done > 'mito_wolb_cov.stats'
```

Now the data is downloaded into R and generate some plots and stats

```R
# load libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)

#first, I have to read the nuclear cov stat and estimate the mean and sd
#then, to add it to 'mito_wolb_cov.stats

nuc_mito_wb_cov <- read.table('data/mito_wolb_cov.stats', header = F) %>% as_tibble()

colnames(nuc_mito_wb_cov) <- c('ID', 'nuc_cov', 'sd_nuc_cov', 'mito_cov', 'wb_cov')
nuc_mito_wb_cov$ID <- str_replace(nuc_mito_wb_cov$ID, '.merged', '')

write_csv(nuc_mito_wb_cov, 'nuc_mit_wb_cov.csv')

# nuclear, mitochondrial and Wb DNA coverage ratio

n_m <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=mito_cov/nuc_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Nuc. to mito. genome coverage ratio", y = "Coverage Ratio")

n_wb <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=wb_cov/nuc_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Nuc. to Wolb. genome coverage ratio", y = "Coverage Ratio")

m_wb <- ggplot(nuc_mito_wb_cov, aes(x=ID, y=mito_cov/wb_cov)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Mito. to Wolb. genome coverage ratio", y = "Coverage Ratio")

ggarrange(n_m, n_wb, m_wb, ncol = 3)
ggsave("cov_ratios.png", height=6, width=15)


# list file names
file_names.window <- list.files(path = "data/",pattern = ".merged.100000_window.cov")

# load data using file names, and make a formatted data frame
setwd("~/R/diro_newgenome/cov_stats/data")

data <- purrr::map_df(file_names.window, function(x) {
  data <- read.delim(x, header = F, sep="\t")
  data$V1 <- str_replace(data$V1, 'dirofilaria_immitis_', '')
  data <- tibble::rowid_to_column(data, "NUM")
  cbind(sample_name = gsub(".merged.100000_window.cov","",x), data)
})
colnames(data) <- c("ID", "NUM", "CHR", "START", "END", 
                    "RAW_COVERAGE", "PROPORTION_COVERAGE")

# remove scaffolds, mitochondrial and wolbachia genome
data_nuc <- dplyr::filter(data, !grepl("scaffold|MtDNA|Wb",CHR))
# data$SEX <- str_detect(data$SCAF,"Trichuris_trichiura_1_")
setwd("~/R/diro_newgenome/cov_stats")

# plot the general cov for each sample
ggplot(data_nuc, aes(NUM, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.5) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("ALL_genomewide_coverage_allsamples.png", height=11.25, width=15)

# Let's see only the chrX to explore the sex of the sample
#Plotting with the chr1 helps to see differences
data_nuc %>%
  filter(., grepl("chrX|chr1",CHR)) %>%
  ggplot(aes(NUM, PROPORTION_COVERAGE/(median(PROPORTION_COVERAGE)), group = ID, col = CHR)) +
  geom_point(size=0.2) +
  labs( x = "Genome position" , y = "Relative coverage per 100kb window") +
  theme_bw() + theme(strip.text.x = element_text(size = 6)) +
  facet_wrap(~ID, scales = "free_y")

ggsave("chrXtochr1_genomewide_coverage_allsamples.png", height=11.25, width=15)
```

