#! /bin/bash

HOST_ACCESSION=GCF_000001735.4
BACTERIA_ACCESSION=GCF_000007805.1

ILLUMINA=$HOME/midden/2021/*-dc3000+arabidopsis-*-brc-downloads

SAMPLES_NAME[0]=48Mt1
SAMPLES_R1[0]=$(ls ${ILLUMINA}/*_48Mt1_*_R1.fastq.gz)
SAMPLES_R2[0]=$(ls ${ILLUMINA}/*_48Mt1_*_R2.fastq.gz)

SAMPLES_NAME[1]=48Mt2
SAMPLES_R1[1]=$(ls ${ILLUMINA}/*_48Mt2_*_R1.fastq.gz)
SAMPLES_R2[1]=$(ls ${ILLUMINA}/*_48Mt2_*_R2.fastq.gz)

SAMPLES_NAME[2]=48Mt3
SAMPLES_R1[2]=$(ls ${ILLUMINA}/*_48Mt3_*_R1.fastq.gz)
SAMPLES_R2[2]=$(ls ${ILLUMINA}/*_48Mt3_*_R2.fastq.gz)

SAMPLES_NAME[3]=48Wt1
SAMPLES_R1[3]=$(ls ${ILLUMINA}/*_48Wt1_*_R1.fastq.gz)
SAMPLES_R2[3]=$(ls ${ILLUMINA}/*_48Wt1_*_R2.fastq.gz)

SAMPLES_NAME[4]=48Wt2
SAMPLES_R1[4]=$(ls ${ILLUMINA}/*_48Wt2_*_R1.fastq.gz)
SAMPLES_R2[4]=$(ls ${ILLUMINA}/*_48Wt2_*_R2.fastq.gz)

SAMPLES_NAME[5]=48Wt3
SAMPLES_R1[5]=$(ls ${ILLUMINA}/*_48Wt3_*_R1.fastq.gz)
SAMPLES_R2[5]=$(ls ${ILLUMINA}/*_48Wt3_*_R2.fastq.gz)

ADDITIONAL_BACTERIA_GENES=local/srna.gtf

# ------------------------------------------------------------------------

if [ -e /programs/docker/bin/docker1 ] ; then
    export HOWTO_DOCKER_CMD=/programs/docker/bin/docker1
fi

# if DONT_USE_STUBS is set to any non-empty string, then don't use
# ./stubs/
#DONT_USE_STUBS=yes

# uncomment to use conda
#CONDA_PREFIX=$(dirname $(dirname $(type -p conda)))
#. "${CONDA_PREFIX}/etc/profile.d/conda.sh"
#conda activate

# Override the default number of threads (nproc --all)
#THREADS=32

