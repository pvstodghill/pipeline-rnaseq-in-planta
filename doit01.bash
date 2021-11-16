#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 1. FASTQC, round 1
# ------------------------------------------------------------------------

echo 1>&2 '# Running FASTQC on raw reads'

rm -rf ${FASTQC1}

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${FASTQC1}/${SAMPLES_NAME[$i]}_R1
    mkdir -p ${FASTQC1}/${SAMPLES_NAME[$i]}_R1
    fastqc -t ${THREADS} \
	   -o ${FASTQC1}/${SAMPLES_NAME[$i]}_R1 \
	   ${INPUTS}/raw_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FASTQC1}/${SAMPLES_NAME[$i]}_R2
	mkdir -p ${FASTQC1}/${SAMPLES_NAME[$i]}_R2
	fastqc -t ${THREADS} \
	       -o ${FASTQC1}/${SAMPLES_NAME[$i]}_R2 \
	       ${INPUTS}/raw_${i}_R2.fastq.gz
    fi
done

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
