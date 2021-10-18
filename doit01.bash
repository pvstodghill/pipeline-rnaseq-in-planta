#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs

# ------------------------------------------------------------------------
# Step 1. FASTQC, round 1
# ------------------------------------------------------------------------

echo 1>&2 '# Running FASTQC on raw reads'

FASTQC1=data/01_fastqc

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
# Step 2. Run FASTP on Illumina reads
# ------------------------------------------------------------------------

echo 1>&2 '# Clean-up Illumina reads...'

FASTP=data/02_fastp

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
# Step 3. FASTQC, round 2
# ------------------------------------------------------------------------

if [ "$RUN_FASTQC_AFTER_FASTP" ] ; then

    echo 1>&2 '# Running FASTQC on trimmed reads'

    FASTQC2=data/03_fastqc

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

fi


# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'

