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