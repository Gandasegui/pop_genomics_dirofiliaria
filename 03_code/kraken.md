# Let's kraken the files to assess the contamination rate

```bash
#loading the module
module load kraken2/2.0.8-beta

#The first step is to build a constum database
WORKING_DIR=/lustre/scratch118/infgen/team333/jg34/POP_Diro
mkdir ${WORKING_DIR}/05_ANALYSIS/KRAKEN
cd KRAKEN

DBNAME='/scratch/c.c2097288/di_project/03_TRIM/KRAKEN/kraken_dog_Di.db'

#First, the taxonomi is downloaded
kraken2-build --download-taxonomy --db $DBNAME

#Now, we will download the ref genome from Canis lupus (dog), upload to hawk and gunzip it
#For adding it to the kraken databse:
kraken2-build --add-to-library GCF_014441545.1_ROS_Cfam_1.0_genomic.fna --no-masking --db $DBNAME
```

--no-masking flag is added becasue dustmasker is not available

so, Masking of Low-complexity Sequences is not possible

```bash
#Now, we also add Di genome
#We have dowloanded it from NCBI databes, since it includes the taxonomy info from NCBI
kraken2-build --add-to-library GCA_001077395.1_ASM107739v1_genomic.fna --no-masking --db $DBNAME

#Now, we build the database
kraken2-build --build --db $DBNAME

# Let's create the while loop to kraken all this things
# The loop:
while read line; do
  kraken2 --db $DBNAME --report $line\.kraken --paired ../$line\/$line\_val_1.fq.gz ../$line\/$line\_val_2.fq.gz;
done < $filename
```

Where filename is a list with all the fastq files used in this work

```bash
#Now, we are going to organize the working directory and remane the samples according to our code
mkdir kraken_report
cp *.kraken kraken_report
cd kraken_report

# Let's recode the samples
mv HW3_1.kraken AUS_QUE_MFP_001-1.kraken
mv HW3_2.kraken AUS_QUE_MFP_001-2.kraken
mv HW4_1.kraken AUS_QUE_MFP_002-1.kraken
mv HW4_2.kraken AUS_QUE_MFP_002-2.kraken
mv HW5_1.kraken AUS_QUE_MFP_003-1.kraken
mv HW5_2.kraken AUS_QUE_MFP_003-2.kraken
mv HW6_1.kraken AUS_QUE_MFP_004-1.kraken
mv HW6_2.kraken AUS_QUE_MFP_004-2.kraken
mv HW6_3.kraken AUS_QUE_MFP_004-3.kraken
mv HW6_4.kraken AUS_QUE_MFP_004-4.kraken
mv HW7_1.kraken AUS_QUE_MFP_005-1.kraken
mv HW7_2.kraken AUS_QUE_MFP_005-2.kraken
mv HW8_1.kraken AUS_QUE_MFP_006-1.kraken
mv HW8_2.kraken AUS_QUE_MFP_006-2.kraken
mv HW9_1.kraken AUS_QUE_MFP_007-1.kraken
mv HW9_2.kraken AUS_QUE_MFP_007-2.kraken
mv HW9_3.kraken AUS_QUE_MFP_007-3.kraken
mv HW12_1.kraken AUS_QUE_MFP_008-1.kraken
mv HW12_2.kraken AUS_QUE_MFP_008-2.kraken
mv HW12_3.kraken AUS_QUE_MFP_008-3.kraken
mv HW13_1.kraken AUS_QUE_MFP_009-1.kraken
mv HW13_2.kraken AUS_QUE_MFP_009-2.kraken
mv HW14_1.kraken AUS_QUE_MFP_010-1.kraken
mv HW14_2.kraken AUS_QUE_MFP_010-2.kraken
mv HW14_3.kraken AUS_QUE_MFP_010-3.kraken
mv HW15_1.kraken AUS_QUE_MFP_011-1.kraken
mv HW15_2.kraken AUS_QUE_MFP_011-2.kraken
mv HW15_3.kraken AUS_QUE_MFP_011-3.kraken
mv HW15_4.kraken AUS_QUE_MFP_011-4.kraken
mv HW16_1.kraken AUS_QUE_MFP_012-1.kraken
mv HW16_2.kraken AUS_QUE_MFP_012-2.kraken
mv GenePool.kraken USA_GEO_ADF_001.kraken
mv Basel_1.kraken	ITL_PAV_ADF_001-1.kraken
mv Basel_2.kraken	ITL_PAV_ADF_001-2.kraken
mv Basel_3.kraken	ITL_PAV_ADF_001-3.kraken
mv DiYaz.kraken USA_MIP_MFS_001.kraken
mv DiJYD.kraken USA_ILL_MFS_001.kraken
mv DiMet.kraken USA_LOU_MFS_001.kraken
mv DiMO.kraken USA_GEO_MFS_001.kraken
mv DiGA2.kraken USA_GEO_MFS_002.kraken
mv JS3.kraken AUS_SYD_AD_001.kraken
mv JS4.kraken AUS_SYD_AD_002.kraken
mv JS5.kraken AUS_SYD_AD_003.kraken
mv JS6.kraken AUS_SYD_AD_004.kraken
mv JS7.kraken AUS_SYD_AD_005.kraken
mv 14B.kraken USA_TEN_MFP_001.kraken
mv 37A.kraken USA_TEX_MFP_001.kraken
mv 39A.kraken USA_GEO_MFP_001.kraken
mv 74A.kraken USA_ARK_MFP_001.kraken
mv 77B.kraken USA_MCH_MFP_001.kraken
mv 20A.kraken USA_GEO_MFP_002.kraken
mv 6B.kraken USA_LOU_MFP_001.kraken

# Now, we multiqc the kraken files
module load multiqc/1.10.1
multiqc *.kraken --title kraken
cd ..
```

Also, we used the 'MiniKraken2_v1_8GB' databse from the kraken website

```bash
while read line; do

kraken2 --db minikraken2_v1_8GB --report $line\.minikraken --paired ../$line\/$line\_val_1.fq.gz ../$line\/$line\_val_2.fq.gz ;

done < $filename

#Let's organized the work enviroment similarly

mkdir minikraken_report
cp *.minikraken minikraken_report
cd minikraken_report

# Let's recode the samples
mv HW3_1.minikraken AUS_QUE_MFP_001-1.minikraken
mv HW3_2.minikraken AUS_QUE_MFP_001-2.minikraken
mv HW4_1.minikraken AUS_QUE_MFP_002-1.minikraken
mv HW4_2.minikraken AUS_QUE_MFP_002-2.minikraken
mv HW5_1.minikraken AUS_QUE_MFP_003-1.minikraken
mv HW5_2.minikraken AUS_QUE_MFP_003-2.minikraken
mv HW6_1.minikraken AUS_QUE_MFP_004-1.minikraken
mv HW6_2.minikraken AUS_QUE_MFP_004-2.minikraken
mv HW6_3.minikraken AUS_QUE_MFP_004-3.minikraken
mv HW6_4.minikraken AUS_QUE_MFP_004-4.minikraken
mv HW7_1.minikraken AUS_QUE_MFP_005-1.minikraken
mv HW7_2.minikraken AUS_QUE_MFP_005-2.minikraken
mv HW8_1.minikraken AUS_QUE_MFP_006-1.minikraken
mv HW8_2.minikraken AUS_QUE_MFP_006-2.minikraken
mv HW9_1.minikraken AUS_QUE_MFP_007-1.minikraken
mv HW9_2.minikraken AUS_QUE_MFP_007-2.minikraken
mv HW9_3.minikraken AUS_QUE_MFP_007-3.minikraken
mv HW12_1.minikraken AUS_QUE_MFP_008-1.minikraken
mv HW12_2.minikraken AUS_QUE_MFP_008-2.minikraken
mv HW12_3.minikraken AUS_QUE_MFP_008-3.minikraken
mv HW13_1.minikraken AUS_QUE_MFP_009-1.minikraken
mv HW13_2.minikraken AUS_QUE_MFP_009-2.minikraken
mv HW14_1.minikraken AUS_QUE_MFP_010-1.minikraken
mv HW14_2.minikraken AUS_QUE_MFP_010-2.minikraken
mv HW14_3.minikraken AUS_QUE_MFP_010-3.minikraken
mv HW15_1.minikraken AUS_QUE_MFP_011-1.minikraken
mv HW15_2.minikraken AUS_QUE_MFP_011-2.minikraken
mv HW15_3.minikraken AUS_QUE_MFP_011-3.minikraken
mv HW15_4.minikraken AUS_QUE_MFP_011-4.minikraken
mv HW16_1.minikraken AUS_QUE_MFP_012-1.minikraken
mv HW16_2.minikraken AUS_QUE_MFP_012-2.minikraken
mv GenePool.minikraken USA_GEO_AD_001.minikraken
mv Basel_1.minikraken	ITL_PAV_ADF_001-1.minikraken
mv Basel_2.minikraken	ITL_PAV_ADF_001-2.minikraken
mv Basel_3.minikraken	ITL_PAV_ADF_001-3.minikraken
mv DiYaz.minikraken USA_MIP_MFS_001.minikraken
mv DiJYD.minikraken USA_ILL_MFS_001.minikraken
mv DiMet.minikraken USA_LOU_MFS_001.minikraken
mv DiMO.minikraken USA_GEO_MFS_001.minikraken
mv DiGA2.minikraken USA_GEO_MFS_002.minikraken
mv JS3.minikraken AUS_SYD_AD_001.minikraken
mv JS4.minikraken AUS_SYD_AD_002.minikraken
mv JS5.minikraken AUS_SYD_AD_003.minikraken
mv JS6.minikraken AUS_SYD_AD_004.minikraken
mv JS7.minikraken AUS_SYD_AD_005.minikraken
mv 14B.minikraken USA_TEN_MFP_001.minikraken
mv 37A.minikraken USA_TEX_MFP_001.minikraken
mv 39A.minikraken USA_GEO_MFP_001.minikraken
mv 74A.minikraken USA_ARK_MFP_001.minikraken
mv 77B.minikraken USA_MCH_MFP_001.minikraken
mv 20A.minikraken USA_GEO_MFP_002.minikraken
mv 6B.minikraken USA_LOU_MFP_001.minikraken

# Now, we multiqc the kraken files
multiqc *.minikraken --title minikraken
```
### The .htlm files can be downloaded and analised using chrome
