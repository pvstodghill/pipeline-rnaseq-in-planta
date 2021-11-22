#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 2. Run FASTP on Illumina reads
# ------------------------------------------------------------------------

echo 1>&2 '# Clean-up Illumina reads...'

rm -rf ${FASTP}
mkdir -p ${FASTP}

# FIXME
## Illumina Stranded
# FASTP_PE_ARGS="--trim_front1 1 --trim_front2 1 --adapter_sequence CTGTCTCTTATACACATCT"
# FASTP_SE_ARGS="--trim_front1 1 --adapter_sequence CTGTCTCTTATACACATCT"
## QIAseq Stranded
FASTP_PE_ARGS=""
FASTP_SE_ARGS=""


for i in $SAMPLES_INDICES ; do
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FASTP}/trimmed_${i}_R2.fastq.gz
	fastp \
	    --thread ${THREADS} \
	    ${FASTP_PE_ARGS} \
	    --json ${FASTP}/${SAMPLES_NAME[$i]}.json \
	    --html ${FASTP}/${SAMPLES_NAME[$i]}.html \
	    --in1 ${INPUTS}/raw_${i}_R1.fastq.gz \
	    --in2 ${INPUTS}/raw_${i}_R2.fastq.gz \
	    --out1 ${FASTP}/trimmed_${i}_R1.fastq.gz \
	    --out2 ${FASTP}/trimmed_${i}_R2.fastq.gz \
	    --unpaired1 ${FASTP}/unpaired_$i.fastq.gz \
	    --unpaired2 ${FASTP}/unpaired_$i.fastq.gz
    else
	echo 1>&2 '##' $i':' ${FASTP}/trimmed_${i}_R1.fastq.gz
	fastp \
	    --thread ${THREADS} \
	    ${FASTP_SE_ARGS} \
	    --json ${FASTP}/${SAMPLES_NAME[$i]}.json \
	    --html ${FASTP}/${SAMPLES_NAME[$i]}.html \
	    --in1 ${INPUTS}/raw_${i}_R1.fastq.gz \
	    --out1 ${FASTP}/trimmed_${i}_R1.fastq.gz
    fi
done

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
