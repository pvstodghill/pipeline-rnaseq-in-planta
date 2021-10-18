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

# ------------------------------------------------------------------------
# Step 4. Make count tables
# ------------------------------------------------------------------------

echo 1>&2 '# Making count tables'

COUNTS=data/04_counts
rm -rf ${COUNTS}
mkdir -p ${COUNTS}

FEATURECOUNTS_ARGS=
FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level

FEATURECOUNTS_ARGS+=" -O" # Assign reads to all their overlapping features

#FEATURECOUNTS_ARGS+=" -M" # all multi-mapping reads reported alignments will be counted
##FEATURECOUNTS_ARGS+=" --fraction" # Assign fractional counts to features
#FEATURECOUNTS_ARGS+=" --primary" # Count primary alignments only !(0x100)

FEATURECOUNTS_ARGS+=" -s 1" # stranded
#FEATURECOUNTS_ARGS+=" -s 2" # reverse-stranded


FEATURECOUNTS_ARGS+=" -p" # fragments (or pairs) will be counted instead of reads (<v2.0.2)
#FEATURECOUNTS_ARGS+=" -p" # input data contains paired-end reads. (>=v2.0.2)
#FEATURECOUNTS_ARGS+=" --countReadPairs" # Count read pairs (fragments) instead of reads

FEATURECOUNTS_ARGS+=" -B" # Only count read pairs that have both ends aligned.
FEATURECOUNTS_ARGS+=" -P" # Check validity of paired-end distance
FEATURECOUNTS_ARGS+=" -C" # Only count concordant reads

FEATURECOUNTS_ARGS+=" -T ${THREADS}"

FEATURECOUNTS_ARGS+=" "
FEATURECOUNTS_ARGS+=" "

fgrep $'\t'gene$'\t' ${INPUTS}/annotation-host.gtf \
      | ./scripts/sanitize-gtf-for-featureCounts \
	    > ${COUNTS}/annotation-host.gtf
fgrep $'\t'gene$'\t' ${INPUTS}/annotation-bacteria.gtf \
      | ./scripts/sanitize-gtf-for-featureCounts \
	    > ${COUNTS}/annotation-bacteria.gtf
cat ${COUNTS}/annotation-host.gtf ${COUNTS}/annotation-bacteria.gtf \
    > ${COUNTS}/annotation-both.gtf


for i in $SAMPLES_INDICES ; do

    echo 1>&2 '##' $i':' ${COUNTS}/both_$i.txt
    featureCounts $FEATURECOUNTS_ARGS \
    		  -a ${COUNTS}/annotation-both.gtf \
    		  -o ${COUNTS}/both_$i.txt \
    		  ${STAR}/aligned_$i.bam

    echo 1>&2 '##' $i':' ${COUNTS}/host_$i.txt
    featureCounts $FEATURECOUNTS_ARGS \
		  -a ${COUNTS}/annotation-host.gtf \
		  -o ${COUNTS}/host_$i.txt \
		  ${STAR}/aligned_$i.bam

    echo 1>&2 '##' $i':' ${COUNTS}/bacteria_$i.txt
    featureCounts $FEATURECOUNTS_ARGS \
		  -a ${COUNTS}/annotation-bacteria.gtf \
		  -o ${COUNTS}/bacteria_$i.txt \
		  ${STAR}/aligned_$i.bam

done


# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'

