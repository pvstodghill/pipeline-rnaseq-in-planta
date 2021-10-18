#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs
FASTP=data/02_fastp

# ------------------------------------------------------------------------
# Step 3. Run Star
# ------------------------------------------------------------------------

echo 1>&2 '# Indexing genomes...'

STAR=data/03_star
rm -rf ${STAR}
mkdir -p ${STAR}

cat ${INPUTS}/genome-host.fna ${INPUTS}/genome-bacteria.fna > ${STAR}/genomes.fna

STAR --runMode genomeGenerate \
     --runThreadN ${THREADS} \
     --genomeDir ${STAR}/genomes \
     --genomeFastaFiles ${STAR}/genomes.fna \
     --genomeSAindexNbases 12 \
     --sjdbGTFfile ${INPUTS}/annotation-host.gtf \
     --sjdbOverhang 49

echo 1>&2 '# Aligning the reads to the genomess...'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${STAR}/aligned_$i.bam

    if [ "$PE" ] ; then

	STAR --runMode alignReads \
	     --genomeDir ${STAR}/genomes \
	     --runThreadN ${THREADS} \
	     --readFilesIn ${FASTP}/trimmed_${i}_R1.fastq.gz ${FASTP}/trimmed_${i}_R2.fastq.gz \
	     --readFilesCommand zcat \
	     --outFilterMultimapNmax 1 \
	     --outFileNamePrefix ${STAR}/${i}. \
	     --outSAMtype BAM Unsorted \
	     --outReadsUnmapped Fastx

    else

	STAR --runMode alignReads \
	     --genomeDir ${STAR}/genomes \
	     --runThreadN ${THREADS} \
	     --readFilesIn ${FASTP}/trimmed_${i}_R1.fastq.gz \
	     --readFilesCommand zcat \
	     --outFilterMultimapNmax 1 \
	     --outFileNamePrefix ${STAR}/${i}. \
	     --outSAMtype BAM Unsorted \
	     --outReadsUnmapped Fastx
    fi


    echo 1>&2 '###' samtools sort aligned_$i.bam
    samtools sort -@ ${THREADS} ${STAR}/${i}.Aligned.out.bam -o ${STAR}/aligned_$i.bam
    rm -f ${STAR}/${i}.Aligned.out.bam

    echo 1>&2 '###' samtools index aligned_$i.bam
    samtools index -@ ${THREADS} ${STAR}/aligned_$i.bam

    echo 1>&2 '###' creating unaligned_${i}_R1.fastq.gz
    mv ${STAR}/${i}.Unmapped.out.mate1 ${STAR}/unaligned_${i}_R1.fastq
    mv ${STAR}/${i}.Unmapped.out.mate2 ${STAR}/unaligned_${i}_R2.fastq
    pigz ${STAR}/unaligned_${i}_R?.fastq

done
