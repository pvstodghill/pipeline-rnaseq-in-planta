#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 3. FASTQC, round 2
# ------------------------------------------------------------------------

echo 1>&2 '# Running FASTQC on trimmed reads'

rm -rf ${FASTQC2}

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${FASTQC2}/${SAMPLES_NAME[$i]}_R1
    mkdir -p ${FASTQC2}/${SAMPLES_NAME[$i]}_R1
    fastqc -t ${THREADS} \
	   -o ${FASTQC2}/${SAMPLES_NAME[$i]}_R1 \
	   ${FASTP}/trimmed_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FASTQC2}/${SAMPLES_NAME[$i]}_R2
	mkdir -p ${FASTQC2}/${SAMPLES_NAME[$i]}_R2
	fastqc -t ${THREADS} \
	       -o ${FASTQC2}/${SAMPLES_NAME[$i]}_R2 \
	       ${FASTP}/trimmed_${i}_R2.fastq.gz
    fi
done



# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
