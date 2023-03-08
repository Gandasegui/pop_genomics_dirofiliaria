# Dirofilaria analysis



REF=dirofilaria_immitis.PRJEB1797.WBPS15.genomic.fa
REF_PREFIX=Di.WBP15
WORKING_DIR=/nfs/users/nfs_s/sd21/lustre118_link/dirofilaria_immitis/POPGEN



DATA_DIR=${WORKING_DIR}/DATA
REF_DIR=${WORKING_DIR}/REF
MAPPING_DIR=${WORKING_DIR}/MAPPING


module load gatk/4.1.4.1
THREADS=20



# CreateSequenceDictionary

[ -f ${REF_DIR}/${REF}.fai ] || samtools faidx ${REF_DIR}/${REF} 

[ -f ${REF_DIR}/${REF%.fa}.dict ] || gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${TMP_DIR}" CreateSequenceDictionary \
     --REFERENCE ${REF_DIR}/${REF} \
     --OUTPUT ${REF_DIR}/${REF%.fa}.dict \
     --spark-runner LOCAL



for FASTQ in `ls -1 ${DATA_DIR}/*_1.fastq.gz`; do

PREFIX=$(basename ${FASTQ} _1.fastq.gz)

if [ -d "${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping" ]; then
	echo -e "\nThere is already a run started with this sample name. Rename and start again\n"
    exit 0
fi

mkdir -p ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp 

PLATFORM_UNIT=$( zcat ${FASTQ} | head -n1 | awk '{for(i=1;i<=NF;i++){if($i~/\:/){a=$i}} print a}' | cut -f1 -d ":" | sed 's/@//g')
DATE=$(date -Iminutes)

# FastqToSam: convert fastq to uBAM
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp -Dsamjdk.compression_level=1" FastqToSam \
     --FASTQ ${DATA_DIR}/${PREFIX}_1.fastq.gz \
     --FASTQ2 ${DATA_DIR}/${PREFIX}_2.fastq.gz \
     --OUTPUT ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.unaligned.bam \
     --SAMPLE_NAME ${PREFIX} \
     --LIBRARY_NAME ${PREFIX}.lib \
     --READ_GROUP_NAME ${PREFIX} \
     --SORT_ORDER queryname \
     --TMP_DIR ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp \
     --spark-runner LOCAL

# MarkIlluminaAdapters: mark adapters in the raw reads
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${PREFIX}_${REF_PREFIX}_mapping/.tmp -Dsamjdk.compression_level=1" MarkIlluminaAdapters \
     --INPUT ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.unaligned.bam \
     --OUTPUT ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.adaptMarked.bam \
     --METRICS ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.adaptMarked.metrics.txt \
     --TMP_DIR ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp \
     --spark-runner LOCAL

# SamToFastq: make interleaved fastq for mapping 
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp -Dsamjdk.compression_level=1" SamToFastq \
     --INPUT ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.adaptMarked.bam \
     --FASTQ ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.interleaved.fastq.gz \
     --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
     --TMP_DIR ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp \
     --spark-runner LOCAL

# mapping
# bwa mem -K 100000000 -v 3 -t ${THREADS} -Y -p ${REF_DIR}/${REF} ${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.interleaved.fastq.gz | samtools view -h -b - > ${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.aligned.bam

minimap2 -t ${THREADS} -a -Y -x sr ${REF_DIR}/${REF} ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.interleaved.fastq.gz | samtools view -h -b - > ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.aligned.bam

# MergeBamAlignment: sort alignment 
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp -Dsamjdk.compression_level=1" MergeBamAlignment \
     --REFERENCE_SEQUENCE ${REF_DIR}/${REF} \
     --UNMAPPED_BAM ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.unaligned.bam \
     --ALIGNED_BAM ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.aligned.bam \
     --OUTPUT ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.alnMerged.bam \
     --CREATE_INDEX false --ADD_MATE_CIGAR true --CLIP_ADAPTERS true --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY BestMapq --ATTRIBUTES_TO_RETAIN XS \
     --spark-runner LOCAL

# MarkDuplicatesSpark: mark duplicates
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/.tmp -Dsamjdk.compression_level=5" MarkDuplicatesSpark \
     --input ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.alnMerged.bam \
     --output ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.deduped.bam \
     --metrics-file ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.deduped.metrics.txt \
     --spark-runner LOCAL

samtools flagstat ${MAPPING_DIR}/${PREFIX}_mapping/${PREFIX}.deduped.bam > ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.deduped.bam.flagstat
samtools stats ${MAPPING_DIR}/${PREFIX}_mapping/${PREFIX}.deduped.bam > ${MAPPING_DIR}/${PREFIX}_${REF_PREFIX}_mapping/${PREFIX}.deduped.bam.stats;

done









# merging datasets
cat ../DATA/lanes_samples.list | while read LANE NAME; do 
     ls $PWD/${NAME}*/*bam | grep -v "unmapped" > ${NAME}.bamlist; 
     done

cat ../DATA/lanes_samples.list | cut -f2 | sort | uniq | while read NAME; do 
     bsub.py --queue long 10 merge "samtools merge -b ${NAME}.bamlist -o ${NAME}.merged.bam";
     done