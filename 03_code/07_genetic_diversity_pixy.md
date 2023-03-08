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

``
