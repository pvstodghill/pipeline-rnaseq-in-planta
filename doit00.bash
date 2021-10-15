#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 0. Set up
# ------------------------------------------------------------------------

rm -rf data

echo 1>&2 '# Initializing data/...'
mkdir -p data/tmp

INPUTS=data/00_inputs
mkdir -p ${INPUTS}

# --------------------------------------------------

(
    cd ${INPUTS}

    echo 1>&2 '## Downloading "datasets" command.'
    curl -s -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets' 
    chmod +x datasets

    echo 1>&2 '## Downloading host and bacteria genomes.'
    
    ./datasets download genome accession --no-progressbar --include-gbff --include-gtf \
	       ${HOST_ACCESSION} ${BACTERIA_ACCESSION}
    unzip -q ncbi_dataset.zip

    echo 1>&2 '## Staging genomes'

    cp ncbi_dataset/data/${HOST_ACCESSION}/GC?_*.fna genome-host.fna
    cp ncbi_dataset/data/${HOST_ACCESSION}/genomic.gff annotation-host.gff
    cp ncbi_dataset/data/${HOST_ACCESSION}/genomic.gtf annotation-host.gtf

    cp ncbi_dataset/data/${BACTERIA_ACCESSION}/GC?_*.fna genome-bacteria.fna
    cp ncbi_dataset/data/${BACTERIA_ACCESSION}/genomic.gff annotation-bacteria.gff
    cp ncbi_dataset/data/${BACTERIA_ACCESSION}/genomic.gtf annotation-bacteria.gtf
)

# --------------------------------------------------

echo 1>&2 '## Making copies of raw reads.'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${INPUTS}/raw_${i}_R1.fastq.gz "<-" "${SAMPLES_R1[$i]}"
    cp "${SAMPLES_R1[$i]}" ${INPUTS}/raw_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${INPUTS}/raw_${i}_R2.fastq.gz "<-" "${SAMPLES_R2[$i]}"
	cp "${SAMPLES_R2[$i]}" ${INPUTS}/raw_${i}_R2.fastq.gz
    fi
done

