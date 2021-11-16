#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 8. Compute some high-level stats
# ------------------------------------------------------------------------

echo 1>&2 '# Generating statistics...'

rm -rf ${STATS}
mkdir -p ${STATS}/temp


cp ${COUNTS}/annotation.gtf ${STATS}/
if [ "${QC_GENES}" ] ; then
    cat "${QC_GENES}" >> ${STATS}/annotation.gtf
fi

# ------------------------------------------------------------------------

init_FEATURECOUNTS_ARGS

FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_biotype"
#FEATURECOUNTS_ARGS+=" -f" # count at feature (exon) level, not the meta-feature (gene) level
FEATURECOUNTS_ARGS+=" -T ${THREADS}"

featureCounts $FEATURECOUNTS_ARGS \
    	      -a ${STATS}/annotation.gtf \
    	      -o ${STATS}/raw_counts.txt \
    	      ${BOWTIE2}/aligned_0.bam \
    	      ${BOWTIE2}/aligned_1.bam \
    	      ${BOWTIE2}/aligned_2.bam \
    	      ${BOWTIE2}/aligned_3.bam \
    	      ${BOWTIE2}/aligned_4.bam \
    	      ${BOWTIE2}/aligned_5.bam


(
    echo -n gene_biotype
    echo -n $'\t'${SAMPLES_NAME[0]}
    echo -n $'\t'${SAMPLES_NAME[1]}
    echo -n $'\t'${SAMPLES_NAME[2]}
    echo -n $'\t'${SAMPLES_NAME[3]}
    echo -n $'\t'${SAMPLES_NAME[4]}
    echo -n $'\t'${SAMPLES_NAME[5]}
    echo ''
    tail -n+3 ${STATS}/raw_counts.txt \
	| cut -f 1,7-
) > ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
