# Dirofilaria immitis populaiton genomics: BUSCOs

### stephen doyle



```bash
cd /nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/BUSCO

conda activate busco_5.4.3
```

## nDi.2.2 / WBP17
```bash
bsub.py --queue long --threads 20 60 busco_di_WBP17_nematoda_odb10 \
    "busco --in dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa --out BUSCO_di_wbp17_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_di_WBP17__eukaryota_odb10 \
    "busco --in dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa --out BUSCO_di_wbp17_genome_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```

## ICBAS_JMDir_1.0
```bash
bsub.py --queue long --threads 20 60 busco_di_ICBAS_JMDir_1.0_nematoda_odb10 \
    "busco --in GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa --out BUSCO_di_ICBAS_JMDir_1.0_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_di_ICBAS_JMDir_1.0_eukaryota_odb10 \
    "busco --in GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa --out BUSCO_di_ICBAS_JMDir_1.0_genome_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"


```


## Brugia malayi

```bash 
bsub.py --queue long --threads 20 60 busco_brugiamalayi_wbp17_nematoda_odb10 \
    "busco --in brugia_malayi.PRJNA10729.WBPS17.genomic.fa --out BUSCO_brugiamalayi_wbp17_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_brugiamalayi_wbp17_eukaryota_odb10 \
    "busco --in brugia_malayi.PRJNA10729.WBPS17.genomic.fa --out BUSCO_brugiamalayi_wbp17_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```

## Onchocerca volvulus

```bash 

ln -s /nfs/users/nfs_s/sd21/lustre_link/REFERENCE_SEQUENCES/onchocerca_volvulus/ONCHO_V4.ref.fa
bsub.py --queue long --threads 20 60 busco_ov_v4_nematoda_odb10 \
    "busco --in ONCHO_V4.ref.fa --out BUSCO_ov_v4_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_ov_v4_eukaryota_odb10 \
    "busco --in ONCHO_V4.ref.fa --out BUSCO_ov_v4_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_ov_v3_nematoda_odb10 \
    "busco --in onchocerca_volvulus.PRJEB513.WBPS17.genomic.fa --out BUSCO_ov_v3_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_ov_v3_eukaryota_odb10 \
    "busco --in onchocerca_volvulus.PRJEB513.WBPS17.genomic.fa --out BUSCO_ov_v3_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"





```





## di ragtag

```bash 
bsub.py --queue long --threads 20 60 busco_di_ragtag_nematoda_odb10 \
    "busco --in ragtag.scaffold.fasta --out BUSCO_di_ragtag_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_di_ragtag_eukaryota_odb10 \
    "busco --in ragtag.scaffold.fasta --out BUSCO_di_ragtag_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```


## dimmitis_WSI_2.2
```bash
ln -s ../dimmitis_WSI_2.2.fa .

bsub.py --queue long --threads 20 60 busco_dimmitis_WSI_2.2_nematoda_odb10 \
    "busco --in dimmitis_WSI_2.2.fa --out BUSCO_dimmitis_WSI_2.2_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"

bsub.py --queue long --threads 20 60 busco_dimmitis_WSI_2.2_nematoda_odb10_augustus \
    "busco --in dimmitis_WSI_2.2.fa --out BUSCO_dimmitis_WSI_2.2_genome_nematoda_odb10_augustus --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r --augustus"

bsub.py --queue long --threads 20 60 busco_dimmitis_WSI_2.2_nematoda_odb10_augustus_long \
    "busco --in dimmitis_WSI_2.2.fa --out BUSCO_dimmitis_WSI_2.2_genome_nematoda_odb10_augustus_long --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r --augustus --long"

bsub.py --queue long --threads 20 60 busco_dimmitis_WSI_2.2_eukaryota_odb10 \
    "busco --in dimmitis_WSI_2.2.fa --out BUSCO_dimmitis_WSI_2.2_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"
```





```
bsub.py --queue long --threads 20 60 busco_di_merged_nematoda_odb10 \
    "busco --in merged.fa --out BUSCO_di_merged_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_ICBAS_plus_recovered_Di_nematoda_odb10     "busco --in ICBAS_plus_recovered_Di.fa --out BUSCO_ICBAS_plus_recovered_Di_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"
```


## Missing BUSCOs

```

ln -s ../../BUSCO/BUSCO_di_wbp17_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv nDi.2.2_full_table.tsv

ln -s ../../BUSCO/BUSCO_di_ICBAS_JMDir_1.0_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv ICBAS_JMDir_1.0_full_table.tsv

ln -s ../../BUSCO/BUSCO_dimmitis_WSI_2.2_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv WSI_2.2_full_table.tsv

ln -s ../GENOME_MERGE/BUSCO_di_merged_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv merged_full_table.tsv

ln -s ../RAGTAG_OV_DI/BUSCO_ICBAS_plus_recovered_Di_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv ICBAS_plus_recovered_full_table.tsv



for i in *_full_table.tsv; do
    echo ${i%_full_table.tsv} > ${i%full_table.tsv}_missing.txt
    cat ${i} | cut -f1 | grep -v "#" | sort | uniq | while read ID; do 
        grep -m1 "${ID}" ${i} | awk '{if($2=="Missing") print "1"; else print "0"}'; 
        done >> ${i%full_table.tsv}_missing.txt; 
    done

echo "ID" > row_ids.txt
cat nDi.2.2_full_table.tsv | cut -f1 | grep -v "#" | sort | uniq >> row_ids.txt



paste row_ids.txt nDi.2.2__missing.txt ICBAS_JMDir_1.0__missing.txt ICBAS_plus_recovered__missing.txt merged__missing.txt WSI_2.2__missing.txt > missing_data.txt

paste row_ids.txt nDi.2.2__missing.txt ICBAS_JMDir_1.0__missing.txt > missing_data.txt2

paste row_ids.txt ICBAS_JMDir_1.0__missing.txt ICBAS_plus_recovered__missing.txt merged__missing.txt WSI_2.2__missing.txt > missing_data.txt

R

library(UpSetR)

data <- read.table("missing_data.txt", header=T)


awk '{print $1,$2,$3,$4,$5}' OFS="\t" Orthogroups.GeneCount.tsv | awk 'NR>1 {for(i=2;i<=NF;i++)if($i>0)$i=1}1' OFS="\t" > orthogroups.data

Orthogroup	PREVIOUS	PUBLISHED	ce.proteins.unique	hc_new_annotation.proteins.unique
OG0000000	0	0	1	0
OG0000001	0	0	1	0
OG0000002	0	0	1	0
OG0000003	0	0	1	0
OG0000004	1	1	0	1