# SCRITP TO CALL AND FILTER THE VARIANTS

### Creating environmental variables and loading modules

```bash
# working dir
WORKING_DIR=/lustre/scratch118/infgen/team333/jg34/POP_Diro
# load gatk
module load gatk/4.1.4.1
#vcftools
module load vcftools/0.1.16-c4
# also need htslib for tabix
module load common-apps/htslib/1.9.229
```

## Step 1. make GVCFs per sample

```bash
#new folder
mkdir ${WORKING_DIR}/04_VARIANTS/GVCFS
cd ${WORKING_DIR}/04_VARIANTS/GVCFS
# create bam list using full path to bams - this allows bams to be anywhere
ls ${WORKING_DIR}/03_MAP/*.merged.bam > ${WORKING_DIR}/04_VARIANTS/bam.list
BAM_LIST=${WORKING_DIR}/04_VARIANTS/bam.list
#indexing the ref
REFERENCE=${WORKING_DIR}/01_REF/dimmitis_WSI_2.2.fa
REF_DIR=${WORKING_DIR}/01_REF
gatk-launch CreateSequenceDictionary -R ${REFERENCE}
# make a sequences list to allow splitting jobs per scaffold/contig
grep ">" ${REFERENCE} | sed -e 's/>//g' > ${WORKING_DIR}/04_VARIANTS/sequences.list
SEQUENCE=${WORKING_DIR}/04_VARIANTS/sequences.list
ulimit -c unlimited
#generate de .dict file from the ref
gatk --java-options "-Xmx8g -Xms4g" CreateSequenceDictionary \
     --REFERENCE ${REFERENCE} \
     --OUTPUT ${REF_DIR}/dimmitis_WSI_2.2.dict \
     --spark-runner LOCAL
#generate, run and submit the jobs
while read BAM; do \
	n=1
	SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.trimmed.bam//g' )
	mkdir ${SAMPLE}_GATK_HC_GVCF
	mkdir ${SAMPLE}_GATK_HC_GVCF/LOGFILES
	echo "gatk GatherVcfsCloud \\" > ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf
	while read SEQUENCE; do
	echo -e "gatk HaplotypeCaller \\
          --input ${BAM} \\
          --output ${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.tmp.gvcf.gz \\
          --reference ${REFERENCE} \\
          --intervals ${SEQUENCE} \\
          --heterozygosity 0.015 \\
          --indel-heterozygosity 0.01 \\
          --annotation DepthPerAlleleBySample --annotation Coverage --annotation ExcessHet --annotation FisherStrand --annotation MappingQualityRankSumTest --annotation StrandOddsRatio --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation DepthPerSampleHC --annotation QualByDepth \\
          --min-base-quality-score 20 --minimum-mapping-quality 30 --standard-min-confidence-threshold-for-calling 30 \\
          --emit-ref-confidence GVCF " > ${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.${SEQUENCE}.tmp.job_${n};
	echo -e "--input ${PWD}/${SAMPLE}_GATK_HC_GVCF/${n}.${SAMPLE}.${SEQUENCE}.tmp.gvcf.gz \\" >> ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf;
	let "n+=1"; done < ${WORKING_DIR}/04_VARIANTS/sequences.list;
	echo -e "--output ${PWD}/${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.gvcf.gz; tabix -p vcf ${PWD}/${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.gvcf.gz" >> ${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf;
	echo -e "rm ${PWD}/${SAMPLE}_GATK_HC_GVCF/*.tmp.* && \\
          mv ${PWD}/${SAMPLE}_GATK_HC_GVCF/*.[oe] ${PWD}/${SAMPLE}_GATK_HC_GVCF/LOGFILES && \\
          cd ${PWD} && \\
          mv ${PWD}/${SAMPLE}_GATK_HC_GVCF ${PWD}/${SAMPLE}_GATK_HC_GVCF_complete" > ${SAMPLE}_GATK_HC_GVCF/run_clean_${SAMPLE};

	chmod a+x ${SAMPLE}_GATK_HC_GVCF/run_*
	# setup job conditions
	JOBS=$( ls -1 ${SAMPLE}_GATK_HC_GVCF/run_hc_* | wc -l )
	ID="U$(date +%s)"
	#submit job array to call variants put scaffold / contig
	bsub -q long -R'span[hosts=1] select[mem>15000] rusage[mem=15000]' -n 6 -M15000 -J GATK_HC_${ID}_[1-${JOBS}]%100 -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_[1-${JOBS}].e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_[1-${JOBS}].o "./${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.*job_\$LSB_JOBINDEX"
	#submit job to gather gvcfs into a single, per sample gvcf
	bsub -q normal -w "done(GATK_HC_${ID}_[1-$JOBS])" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J GATK_HC_${ID}_gather_gvcfs -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_gather_gvcfs.e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_gather_gvcfs.o "./${SAMPLE}_GATK_HC_GVCF/run_gather_${SAMPLE}_gvcf"
	# clean up
	bsub -q normal -w "done(GATK_HC_${ID}_gather_gvcfs)" -R'span[hosts=1] select[mem>500] rusage[mem=500]' -n 1 -M500 -J GATK_HC_${ID}_clean -e ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_clean.e -o ${SAMPLE}_GATK_HC_GVCF/GATK_HC_${ID}_clean.o "./${SAMPLE}_GATK_HC_GVCF/run_clean_${SAMPLE}"
	sleep 1
done < ${BAM_LIST}
```

## Step 2. Gather the GVCFs to generate a merged GVCF

```bash
# make a new directory for the merged GVCFS
mkdir ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED
cd ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED
# make a list of GVCFs to be merged
ls -1 ${WORKING_DIR}/04_VARIANTS/GVCFS/*complete/*gz > ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED/gvcf.list
GVCF_LIST=${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED/gvcf.list
REFERENCE=${WORKING_DIR}/01_REF/dimmitis_WSI_2.2.fa
# setup the run files
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


## Step 3. Split merged GVCF into individual sequences, and then genotype to generate a VCF

```bash
# split each chromosome up into separate jobs, and run genotyping on each individually.
n=1
while read SEQUENCE; do
     echo -e "gatk GenotypeGVCFs \
     -R ${REFERENCE} \
     -V ${SEQUENCE}.cohort.g.vcf.gz \
     --intervals ${SEQUENCE} \
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
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -n 4 -M10000 -J GATK_HC_GENOTYPE_${ID}_[1-$JOBS] -e LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].e -o LOGFILES/GATK_HC_GENOTYPE_${ID}_[1-$JOBS].o "./run_hc_*\$LSB_JOBINDEX"
```

## Step 4. Bring the files together

```bash
# make list of vcfs
ls -1 *.cohort.vcf.gz | sort -n > vcf_files.list
# merge them
vcf-concat --files vcf_files.list > dirofilaria_immitis.cohort.vcf;
     bgzip dirofilaria_immitis.cohort.vcf;
     tabix -p vcf dirofilaria_immitis.cohort.vcf.gz
# clean up
rm run*
rm ^[0-9]*
rm *.g.vcf.gz*
```

## Querying SNP and INDEL QC profiles to determine thresholds for filters
Adapted from https://evodify.com/gatk-in-non-model-organism/

```bash
mkdir ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED/FILTER
cd ${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED/FILTER

# set reference, vcf, and mitochondrial and Wb contig
VCF=${WORKING_DIR}/04_VARIANTS/GATK_HC_MERGED/FILTER/dirofilaria_immitis.cohort.vcf.gz
MIT_CONTIG=dirofilaria_immitis_chrMtDNA
WB_CONTIG=dirofilaria_immitis_chrWb

# select nuclear SNPs
bsub.py 1 select_nuclearSNPs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearSNPs.vcf"

# select nuclear INDELs
bsub.py 1 select_nuclearINDELs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--exclude-intervals ${MIT_CONTIG} \
--exclude-intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.nuclearINDELs.vcf"

# select mitochondrial SNPs
bsub.py 1 select_mitoSNPs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoSNPs.vcf"

# select mitochondrial INDELs
bsub.py 1 select_mitoINDELs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${MIT_CONTIG} \
--output ${VCF%.vcf.gz}.mitoINDELs.vcf"

# select WB SNPs
bsub.py 1 select_WbSNPs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include SNP \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbSNPs.vcf"

# select WB INDELs
bsub.py 1 select_WbINDELs "gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF} \
--select-type-to-include INDEL \
--intervals ${WB_CONTIG} \
--output ${VCF%.vcf.gz}.WbINDELs.vcf"

# make a table of nuclear SNP data
bsub.py 1 select_nuclearSNPs_table "gatk VariantsToTable \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.nuclearSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearSNPs.table"

# make a table of nuclear INDEL data data
bsub.py 1 select_nuclearINDELs_table "gatk VariantsToTable \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.nuclearINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_nuclearINDELs.table"

# make a table of mito SNP data
bsub.py --done "select_mitoSNPs" 1 select_mitoSNPs_table "gatk VariantsToTable \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.mitoSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoSNPs.table"

# make a table of mito INDEL data data
bsub.py --done "select_mitoINDELs"  1 select_mitoINDELs_table "gatk VariantsToTable \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.mitoINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_mitoINDELs.table"

# make a table of Wb SNP data
bsub.py --done "select_WbSNPs" 1 select_mitoSNPs_table "gatk VariantsToTable \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.WbSNPs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbSNPs.table"

# make a table of Wb INDEL data data
bsub.py --done "select_WbINDELs"  1 select_mitoINDELs_table "gatk VariantsToTable \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.WbINDELs.vcf \
--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
--output GVCFall_WbINDELs.table"
```
### Make some density plots of the data and get quantiles in R

```R
# load libraries
library(patchwork)
require(data.table)
library(tidyverse)
library(gridExtra)

VCF_nuclear_snps <- fread('data/GVCFall_nuclearSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_nuclear_snps <- sample_frac(VCF_nuclear_snps, 0.2)
VCF_nuclear_indels <- fread('data/GVCFall_nuclearINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_nuclear_indels <- sample_frac(VCF_nuclear_indels, 0.2)
dim(VCF_nuclear_snps)
dim(VCF_nuclear_indels)
VCF_nuclear <- rbind(VCF_nuclear_snps, VCF_nuclear_indels)
VCF_nuclear$Variant <- factor(c(rep("SNPs", dim(VCF_nuclear_snps)[1]), rep("Indels", dim(VCF_nuclear_indels)[1])))

VCF_mito_snps <- fread('data/GVCFall_mitoSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_mito_indels <- fread('data/GVCFall_mitoINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
dim(VCF_mito_snps)
dim(VCF_mito_indels)
VCF_mito <- rbind(VCF_mito_snps, VCF_mito_indels)
VCF_mito$Variant <- factor(c(rep("SNPs", dim(VCF_mito_snps)[1]), rep("Indels", dim(VCF_mito_indels)[1])))

VCF_wb_snps <- fread('data/GVCFall_WbSNPs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
VCF_wb_indels <- fread('data/GVCFall_WbINDELs.table', header = TRUE, fill=TRUE, na.strings=c("","NA"), sep = "\t")
dim(VCF_wb_snps)
dim(VCF_wb_indels)
VCF_wb <- rbind(VCF_wb_snps, VCF_wb_indels)
VCF_wb$Variant <- factor(c(rep("SNPs", dim(VCF_wb_snps)[1]), rep("Indels", dim(VCF_wb_indels)[1])))


snps <- '#A9E2E4'
indels <- '#F4CCCA'

fun_variant_summaries <- function(data, title){
  # gatk hardfilter: SNP & INDEL QUAL < 0
  QUAL_quant <- quantile(data$QUAL, c(.01,.99), na.rm=T)
  
  QUAL <-
    ggplot(data, aes(x=log10(QUAL), fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=0, size=0.7, col="red") +
    geom_vline(xintercept=c(log10(QUAL_quant[2]), log10(QUAL_quant[3])), size=0.7, col="blue") +
    #xlim(0,10000) +
    theme_bw() +
    labs(title=paste0(title,": QUAL"))
  
  
  # DP doesnt have a hardfilter
  DP_quant <- quantile(data$DP, c(.01,.99), na.rm=T)
  
  DP <-
    ggplot(data, aes(x=log10(DP), fill=Variant)) +
    geom_density(alpha=0.3) +
    geom_vline(xintercept=log10(DP_quant), col="blue") +
    theme_bw() +
    labs(title=paste0(title,": DP"))
  
  # gatk hardfilter: SNP & INDEL QD < 2
  QD_quant <- quantile(data$QD, c(.01,.99), na.rm=T)
  
  QD <-
    ggplot(data, aes(x=QD, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=2, size=0.7, col="red") +
    geom_vline(xintercept=QD_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": QD"))
  
  # gatk hardfilter: SNP FS > 60, INDEL FS > 200
  FS_quant <- quantile(data$FS, c(.01,.99), na.rm=T)
  
  FS <-
    ggplot(data, aes(x=log10(FS), fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(log10(60), log10(200)), size=0.7, col="red") +
    geom_vline(xintercept=log10(FS_quant), size=0.7, col="blue") +
    #xlim(0,250) +
    theme_bw() +
    labs(title=paste0(title,": FS"))
  
  # gatk hardfilter: SNP & INDEL MQ < 30
  MQ_quant <- quantile(data$MQ, c(.01,.99), na.rm=T)
  
  MQ <-
    ggplot(data, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
    geom_vline(xintercept=40, size=0.7, col="red") +
    geom_vline(xintercept=MQ_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": MQ"))
  
  # gatk hardfilter: SNP MQRankSum < -20
  MQRankSum_quant <- quantile(data$MQRankSum, c(.01,.99), na.rm=T)
  
  MQRankSum <-
    ggplot(data, aes(x=log10(MQRankSum), fill=Variant)) + geom_density(alpha=.3) +
    geom_vline(xintercept=log10(-20), size=0.7, col="red") +
    geom_vline(xintercept=log10(MQRankSum_quant), size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": MQRankSum"))
  
  
  # gatk hardfilter: SNP SOR < 4 , INDEL SOR > 10
  SOR_quant <- quantile(data$SOR, c(.01, .99), na.rm=T)
  
  SOR <-
    ggplot(data, aes(x=SOR, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(4, 10), size=1, colour = c(snps,indels)) +
    geom_vline(xintercept=SOR_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": SOR"))
  
  # gatk hardfilter: SNP ReadPosRankSum <-10 , INDEL ReadPosRankSum < -20
  ReadPosRankSum_quant <- quantile(data$ReadPosRankSum, c(.01,.99), na.rm=T)
  
  ReadPosRankSum <-
    ggplot(data, aes(x=ReadPosRankSum, fill=Variant)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(snps,snps,indels,indels)) +
    xlim(-10, 10) +
    geom_vline(xintercept=ReadPosRankSum_quant, size=0.7, col="blue") +
    theme_bw() +
    labs(title=paste0(title,": ReadPosRankSum"))
  
  
  plot <- QUAL + DP + QD + FS + MQ + MQRankSum + SOR + ReadPosRankSum + plot_layout(ncol=2)
  
  print(plot)
  
  ggsave(paste0("plot_",title,"_variant_summaries.png"), height=20, width=15, type="cairo")
  
  
  # generate a table of quantiles for each variant feature
  QUAL_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(QUAL, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  QUAL_quant$name <- "QUAL"
  DP_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(DP, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  DP_quant$name <- "DP"
  QD_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(QD, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  QD_quant$name <- "QD"
  FS_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(FS, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  FS_quant$name <- "FS"
  MQ_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(MQ, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  MQ_quant$name <- "MQ"
  MQRankSum_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(MQRankSum, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  MQRankSum_quant$name <- "MQRankSum"
  SOR_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(SOR, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  SOR_quant$name <- "SOR"
  ReadPosRankSum_quant <- data %>% group_by(Variant) %>% summarise(quants = list(quantile(ReadPosRankSum, probs = c(0.01,0.05,0.95,0.99),na.rm=T))) %>% unnest_wider(quants)
  ReadPosRankSum_quant$name <- "ReadPosRankSum"
  
  quantiles <- bind_rows(QUAL_quant,DP_quant, QD_quant, FS_quant, MQ_quant, MQRankSum_quant, SOR_quant, ReadPosRankSum_quant)
  quantiles$name <- c("QUAL_Indels","QUAL_SNPs","DP_indels","DP_SNPs", "QD_indels","QD_SNPs", "FS_indels","FS_SNPs", "MQ_indels","MQ_SNPs", "MQRankSum_indels","MQRankSum_SNPs", "SOR_indels","SOR_SNPs","ReadPosRankSum_indels","ReadPosRankSum_SNPs")
  
  png(paste0("table_",title,"_variant_quantiles.png"), width=1000,height=500,bg = "white")
  print(quantiles)
  grid.table(quantiles)
  dev.off()
  
}

# run nuclear variants
fun_variant_summaries(VCF_nuclear,"nuclear")
# run mitochondrial variants
fun_variant_summaries(VCF_mito,"mitochondrial")
# run wolbachia variants
fun_variant_summaries(VCF_wb,"wolbachia")
```
### Attending to the quantiles, thresholds for specific parameters are stablished

```bash
#Nuclear
bsub.py 1 filter_nuclearSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.nuclearSNPs.vcf \
--filter-expression 'QUAL < 36 || DP < 79 || DP > 4154 || MQ < 40.00 || SOR > 5.600 || QD < 0.28 || FS > 13.000 || MQRankSum < -4.800 || ReadPosRankSum < -3.200 || ReadPosRankSum > 2.300' \
--filter-name "SNP_filtered" \
--output dirofilaria_immitis.cohort.nuclearSNPs.filtered.vcf"

bsub.py 1 filter_nuclearINDELs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.nuclearINDELs.vcf \
--filter-expression 'QUAL < 43 || DP < 43 || DP > 3599 || MQ < 44.00 || SOR > 5.800 || QD < 0.73 || FS > 8.000 || MQRankSum < -5.600 || ReadPosRankSum < -4.900 || ReadPosRankSum > 2.400' \
--filter-name "INDEL_filtered" \
--output dirofilaria_immitis.cohort.nuclearINDELs.filtered.vcf"

#Mitochondrial
bsub.py 1 filter_mitoSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.mitoSNPs.vcf \
--filter-expression ' QUAL < 48 || DP < 15003 || DP > 83794 || MQ < 52.00 || SOR > 6.000 || QD < 0.2 || FS > 119.0 || MQRankSum < -13.0 || ReadPosRankSum < -7.8 || ReadPosRankSum > 6.3 ' \
--filter-name "SNP_filtered" \
--output dirofilaria_immitis.cohort.mitoSNPs.filtered.vcf"

bsub.py 1 filter_mitoINDELs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.mitoINDELs.vcf \
--filter-expression 'QUAL < 40 || DP < 14424 || DP > 82781 || MQ < 50.00 || SOR > 2.2000 || QD < 0.7 || FS > 40.0 || ReadPosRankSum < -11.6 || ReadPosRankSum > 8.4' \
--filter-name "INDEL_filtered" \
--output dirofilaria_immitis.cohort.mitoINDELs.filtered.vcf"

#Wolbachia
bsub.py 1 filter_WbSNPs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.WbSNPs.vcf \
--filter-expression ' QUAL < 30 || DP < 907 || DP > 2984 || MQ < 34.00 || SOR > 6.100 || QD < 1.0 || FS > 72.0 || MQRankSum < -8.1 || ReadPosRankSum < -4.2 || ReadPosRankSum > 3.6 ' \
--filter-name "SNP_filtered" \
--output dirofilaria_immitis.cohort.WbSNPs.filtered.vcf"

bsub.py 1 filter_WbINDELs "gatk VariantFiltration \
--reference ${REFERENCE} \
--variant dirofilaria_immitis.cohort.WbINDELs.vcf \
--filter-expression 'QUAL < 47 || DP < 1047 || DP > 3310 || MQ < 46.00 || SOR > 6.000 || QD < 1.5 || FS > 42.0 || ReadPosRankSum < -6.5 || ReadPosRankSum > 2.0' \
--filter-name "INDEL_filtered" \
--output dirofilaria_immitis.cohort.WbINDELs.filtered.vcf"

# once done, count the filtered sites - funny use of "|" allows direct markdown table format
echo -e "| Filtered_VCF | Variants_PASS | Variants_FILTERED |\n| -- | -- | -- | " > filter.stats

for i in *filtered.vcf; do
     name=${i}; pass=$( grep -E 'PASS' ${i} | wc -l ); filter=$( grep -E 'filter' ${i} | wc -l );
     echo -e "| ${name} | ${pass} | ${filter} |" >> filter.stats
done
```
This is the summary of the filtered variants

| Filtered_VCF | Variants_PASS | Variants_FILTERED |
| -- | -- | -- | 
| dirofilaria_immitis.cohort.mitoINDELs.filtered.vcf | 38 | 29 |
| dirofilaria_immitis.cohort.mitoSNPs.filtered.vcf | 88 | 21 |
| dirofilaria_immitis.cohort.nuclearINDELs.filtered.vcf | 492119 | 54848 |
| dirofilaria_immitis.cohort.nuclearSNPs.filtered.vcf | 376359 | 44976 |
| dirofilaria_immitis.cohort.WbINDELs.filtered.vcf | 4262 | 294 |
| dirofilaria_immitis.cohort.WbSNPs.filtered.vcf | 656 | 46 |

```bash

#Now merging all the files
bsub.py 1 merge_nuclear_variants "gatk MergeVcfs \
--INPUT dirofilaria_immitis.cohort.nuclearSNPs.filtered.vcf \
--INPUT dirofilaria_immitis.cohort.nuclearINDELs.filtered.vcf \
--OUTPUT dirofilaria_immitis.cohort.nuclearALL.filtered.vcf"

bsub.py 1 merge_mito_variants "gatk MergeVcfs \
--INPUT dirofilaria_immitis.cohort.mitoSNPs.filtered.vcf \
--INPUT dirofilaria_immitis.cohort.mitoINDELs.filtered.vcf \
--OUTPUT dirofilaria_immitis.cohort.mitoALL.filtered.vcf"

bsub.py 1 merge_mito_variants "gatk MergeVcfs \
--INPUT dirofilaria_immitis.cohort.WbSNPs.filtered.vcf \
--INPUT dirofilaria_immitis.cohort.WbINDELs.filtered.vcf \
--OUTPUT dirofilaria_immitis.cohort.WbALL.filtered.vcf"
```

### Filter genotypes based on x3 depth per genotype

```bash
#Nuclear 
bsub.py 1 filter_nuclear_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearALL.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.nuclearALL.DPfiltered.vcf"
bsub.py --done "filter_nuclear_GT" 1 filter_nuclear_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.nuclearALL.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.nuclearALL.DPfilterNoCall.vcf"

#Mito
bsub.py 1 filter_mito_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoALL.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.mitoALL.DPfiltered.vcf"
bsub.py --done "filter_mito_GT" 1 filter_mito_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.mitoALL.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.mitoALL.DPfilterNoCall.vcf"

#wolbachia
bsub.py 1 filter_Wb_GT \
"gatk VariantFiltration \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbALL.filtered.vcf \
--genotype-filter-expression ' DP < 3 '  \
--genotype-filter-name "DP_lt3" \
--output ${VCF%.vcf.gz}.WbALL.DPfiltered.vcf"
bsub.py --done "filter_Wb_GT" 1 filter_mito_GT2 \
"gatk SelectVariants \
--reference ${REFERENCE} \
--variant ${VCF%.vcf.gz}.WbALL.DPfiltered.vcf \
--set-filtered-gt-to-nocall \
--output ${VCF%.vcf.gz}.WbALL.DPfilterNoCall.vcf"
```

### Now we apply a set of standard filters for population genomics

```bash
#Nuclear variants
vcftools \
--vcf ${VCF%.vcf.gz}.nuclearALL.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--hwe 1e-06 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.nuclear_variants.final
#> After filtering, kept 32 out of 32 Individuals
#> After filtering, kept 529424 out of a possible 968290 Sites
#--- nuclear SNPs
vcftools --vcf dirofilaria_immitis.cohort.nuclear_variants.final.recode.vcf --remove-indels
#> After filtering, kept 351657 out of a possible 529424 Sites
#--- nuclear  INDELs
vcftools --vcf dirofilaria_immitis.cohort.nuclear_variants.final.recode.vcf --keep-only-indels
#> After filtering, kept 177767 out of a possible 529424 Sites

#Mitochondrial variants
vcftools \
--vcf ${VCF%.vcf.gz}.mitoALL.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.mito_variants.final
#> After filtering, kept 32 out of 32 Individuals
#> After filtering, kept 50 out of a possible 164 Sites
#--- mito SNPs
vcftools --vcf dirofilaria_immitis.cohort.mito_variants.final.recode.vcf --remove-indels
#> After filtering, kept 40 out of a possible 50 Sites
#--- mito INDELs
vcftools --vcf dirofilaria_immitis.cohort.mito_variants.final.recode.vcf --keep-only-indels
#> After filtering, kept 10 out of a possible 50 Sites

#wolbachia variants
vcftools \
--vcf ${VCF%.vcf.gz}.WbALL.DPfilterNoCall.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.02 \
--recode \
--recode-INFO-all \
--out ${VCF%.vcf.gz}.Wb_variants.final
#> After filtering, kept 32 out of 32 Individuals
#> After filtering, kept 974 out of a possible 5246 Sites
#--- Wb SNPs
vcftools --vcf dirofilaria_immitis.cohort.Wb_variants.final.recode.vcf --remove-indels
#> After filtering, kept 561 out of a possible 974 Sites
#--- Wb INDELs
vcftools --vcf dirofilaria_immitis.cohort.Wb_variants.final.recode.vcf --keep-only-indels
#> After filtering, kept 413 out of a possible 974 Sites
```
### Now, we are filtering by missingness

```bash
#determine missingness per individual
vcftools --vcf ${VCF%.vcf.gz}.nuclear_variants.final.recode.vcf --out nuclear --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.mito_variants.final.recode.vcf --out mito --missing-indv
vcftools --vcf ${VCF%.vcf.gz}.Wb_variants.final.recode.vcf --out wb --missing-indv
```

### Let's check the missingess in R

```R
data_nuclear <- read.delim("data/nuclear.imiss", header=T)
data_mito <- read.delim("data/mito.imiss", header=T)
data_wb <- read.delim("data/wb.imiss", header=T)
#crearting the function - per sample
fun_plot_missingness <- function(data,title) {
  
  data <- data %>% separate(INDV, c("country","population","sampletype","sampleID"))
  count <- data[1,5]
  
  plot <- ggplot(data, aes(population, 1-F_MISS)) +
    geom_boxplot(color = 'brown') +
    geom_point(size = 1, color = 'brown4') +
    theme_bw() +
    labs(x="Region", y="Proportion of total variants present (1-missingness)")+
    ggtitle(title)
  print(plot)
  ggsave(paste0("plot_missingness_figure",title,".png"))
}
# plotting for each dataset
fun_plot_missingness(data_nuclear, "nuclear_variants")
fun_plot_missingness(data_mito,"mitochondrial_variants")
fun_plot_missingness(data_wb, "wb_variants")

#And now per sample type
fun_plot_missingness_sampletype <- function(data,title) {
  
  data <- data %>% separate(INDV, c("country","population","sampletype","sampleID"))
  data$sampletype <- gsub('ADF', 'AD', data$sampletype)
  count <- data[1,5]
  
  plot <- ggplot(data, aes(sampletype, 1-F_MISS)) +
    geom_boxplot(color = 'royalblue') +
    geom_point(size = 1, color = 'royalblue') +
    theme_bw() +
    labs(x="Country", y="Proportion of total variants present (1-missingness)", title=paste0("Variants per sample: ",title, " (n = ", count,")"))
  print(plot)
  ggsave(paste0("plot_missingness_sampletype_",title,".png"))
}
# plotting for each dataset
fun_plot_missingness_sampletype(data_nuclear, "nuclear_variants")
fun_plot_missingness_sampletype(data_mito,"mitochondrial_variants")
fun_plot_missingness_sampletype(data_wb, "wb_variants")
```

It is variable and as expected, some of the sample from Quesland provided high missingess

So now we will generate different sample list for each database and ten evaluate max missiness

Also, let's remove the Chinise sample

```bash
# For nuclear (n=19) - nuclear_samplelist.keep
AUS_NSW_AD_001
AUS_NSW_AD_002
AUS_NSW_AD_003
AUS_NSW_AD_004
AUS_NSW_AD_005
ITL_PAV_ADF_001
USA_ARK_MFP_001
USA_GEO_ADF_001
USA_GEO_MFP_001
USA_GEO_MFP_002
USA_GEO_MFS_001
USA_GEO_MFS_002
USA_ILL_MFS_001
USA_LOU_MFP_001
USA_LOU_MFS_001
USA_MCH_MFP_001
USA_MIP_MFS_001
USA_TEN_MFP_001
USA_TEX_MFP_001
# For mithochondiral (n=30) - mito_samplelist.keep
AUS_NSW_AD_001
AUS_NSW_AD_002
AUS_NSW_AD_003
AUS_NSW_AD_004
AUS_NSW_AD_005
AUS_QUE_MFP_001
AUS_QUE_MFP_002
AUS_QUE_MFP_003
AUS_QUE_MFP_004
AUS_QUE_MFP_005
AUS_QUE_MFP_006
AUS_QUE_MFP_008
AUS_QUE_MFP_009
AUS_QUE_MFP_010
AUS_QUE_MFP_011
ITL_PAV_ADF_001
USA_ARK_MFP_001
USA_GEO_ADF_001
USA_GEO_MFP_001
USA_GEO_MFP_002
USA_GEO_MFS_001
USA_GEO_MFS_002
USA_ILL_MFS_001
USA_LOU_MFP_001
USA_LOU_MFS_001
USA_MCH_MFP_001
USA_MIP_MFS_001
USA_TEN_MFP_001
USA_TEX_MFP_001
# For wb (n=16) - wb_samplelist.keep
AUS_NSW_AD_001
AUS_NSW_AD_002
AUS_NSW_AD_003
AUS_NSW_AD_004
AUS_NSW_AD_005
ITL_PAV_ADF_001
USA_GEO_ADF_001
USA_GEO_MFP_001
USA_GEO_MFP_002
USA_GEO_MFS_001
USA_GEO_MFS_002
USA_ILL_MFS_001
USA_LOU_MFP_001
USA_LOU_MFS_001
USA_MCH_MFP_001
USA_MIP_MFS_001
```
### Once we selected those samples with low missiness, let's check different thresholds for each dataset

```bash
# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.nuclear_variants.final.recode.vcf --keep nuclear_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 19 out of 32 Individuals
#After filtering, kept 488446 out of a possible 529424 Sites

# max-missing = 0.8
#After filtering, kept 19 out of 32 Individuals
#After filtering, kept 422935 out of a possible 529424 Sites

# max-missing = 0.9
#After filtering, kept 19 out of 32 Individuals
#After filtering, kept 248619 out of a possible 529424 Sites

# max-missing = 1
#After filtering, kept 19 out of 32 Individuals
#After filtering, kept 106209 out of a possible 529424 Sites

# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.mito_variants.final.recode.vcf --keep mito_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 30 out of 32 Individuals
#After filtering, kept 49 out of a possible 50 Sites

# max-missing = 0.8
#After filtering, kept 30 out of 32 Individuals
#After filtering, kept 49 out of a possible 50 Sites

# max-missing = 0.9
#After filtering, kept 30 out of 32 Individuals
#After filtering, kept 46 out of a possible 50 Sites

# max-missing = 1
#After filtering, kept 30 out of 32 Individuals
#After filtering, kept 39 out of a possible 50 Sites

# For nuclear variants
for i in 0.7 0.8 0.9 1; do
     vcftools --vcf ${VCF%.vcf.gz}.Wb_variants.final.recode.vcf --keep wb_samplelist.keep --max-missing ${i} ;
done

# max-missing = 0.7
#After filtering, kept 17 out of 32 Individuals
#After filtering, kept 447 out of a possible 974 Sites

# max-missing = 0.8
#After filtering, kept 17 out of 32 Individuals
#After filtering, kept 73 out of a possible 50 Sites

# max-missing = 0.9
#After filtering, kept 17 out of 32 Individuals
#After filtering, kept 8 out of a possible 50 Sites

# max-missing = 1
#After filtering, kept 17 out of 32 Individuals
#After filtering, kept 0 out of a possible 50 Sites
```
Selecting a max missingness of 0.8 for nuclear and mito, and 0.7 for Wb is sensible

```bash
mkdir ../../FINAL_SETS/
# For nuclear
vcftools --vcf ${VCF%.vcf.gz}.nuclear_variants.final.recode.vcf \
     --keep nuclear_samplelist.keep \
     --max-missing 0.8 \
     --recode --recode-INFO-all \
     --out ../../FINAL_SETS/nuclear_samples3x_missing0.8
# For mito
vcftools --vcf ${VCF%.vcf.gz}.mito_variants.final.recode.vcf \
     --keep mito_samplelist.keep \
     --max-missing 0.8 \
     --recode --recode-INFO-all \
     --out .../../FINAL_SETS/mito_samples3x_missing0.8
# For wb
vcftools --vcf ${VCF%.vcf.gz}.Wb_variants.final.recode.vcf \
     --keep wb_samplelist.keep \
     --max-missing 0.7 \
     --recode --recode-INFO-all \
     --out ../../FINAL_SETS/wb_samples3x_missing0.7
```
### Also, we will slectct only the variants in the chr 1 to chr4, avoising the chrX and the scaffolds

```bash
vcftools --vcf ../FINAL_SETS/nuclear_samples3x_missing0.8.recode.vcf \
--chr dirofilaria_immitis_chr1 \
--chr dirofilaria_immitis_chr2 \
--chr dirofilaria_immitis_chr3 \
--chr dirofilaria_immitis_chr4 \
--recode --out ../FINAL_SETS/nuclear_samples3x_missing0.8.chr1to4
#After filtering, kept 19 out of 19 Individuals
#After filtering, kept 312768 out of a possible 422935 Sites

vcftools --vcf nuclear_samples3x_missing0.8.chr1to4.recode.vcf --remove-indels
#> After filtering, kept 216455 out of a possible 312768 Sites
vcftools --vcf nuclear_samples3x_missing0.8.chr1to4.recode.vcf --keep-only-indels
#> After filtering, kept 96313 out of a possible 312768 Sites
```

### Let's do also a file with no indels for PCA using variant freq

```bash
vcftools --vcf nuclear_samples3x_missing0.8.chr1to4.recode.vcf \
--remove-indels \
--recode --out nuclear_samples3x_missing0.8.chr1to4.NOindels
#> After filtering, kept 216455 out of a possible 312768 Sites
```
