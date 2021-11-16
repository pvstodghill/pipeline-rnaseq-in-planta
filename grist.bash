# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'

#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs
FASTP=data/02_fastp

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'

#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs
FASTP=data/02_fastp
STAR=data/03_star
COUNTS=data/04_counts

# ------------------------------------------------------------------------
# Step 5. Run DESeq2
# ------------------------------------------------------------------------

echo 1>&2 '# Running DESeq2...'

CHANGE_CUTOFF=2.0
TAG=_48Wt-48Mt

DESEQ2=data/05_deseq2
rm -rf ${DESEQ2}
mkdir -p ${DESEQ2}/temp

cp ${COUNTS}/annotation-bacteria.gtf ${DESEQ2}/temp/regions.gtf
./scripts/make-counts-table-from-featurecounts \
    ${DESEQ2}/temp/regions.gtf \
    48Mt1:${COUNTS}/bacteria_0.txt \
    48Mt2:${COUNTS}/bacteria_1.txt \
    48Mt3:${COUNTS}/bacteria_2.txt \
    48Wt1:${COUNTS}/bacteria_3.txt \
    48Wt2:${COUNTS}/bacteria_4.txt \
    48Wt3:${COUNTS}/bacteria_5.txt \
    > ${DESEQ2}/temp/counts.txt

# ------------------------------------------------------------------------

./scripts/prep-deseq2 -x -s ./scripts \
		      -F parametric \
		      -d ${DESEQ2}/temp \
		      -t ${TAG} \
		      -c ${DESEQ2}/temp/counts.txt \
 		     48Wt:48Wt1 48Wt:48Wt2 48Wt:48Wt3 \
 		     48Mt:48Mt1 48Mt:48Mt2 48Mt:48Mt3

# ------------------------------------------------------------------------

Rscript ./scripts/run-deseq2 \
	${DESEQ2}/temp/params${TAG}.R

# ------------------------------------------------------------------------

cat ${DESEQ2}/temp/output-extended${TAG}.txt \
    | ./scripts/deseq-output2results > ${DESEQ2}/results${TAG}.txt

cat ${DESEQ2}/temp/output${TAG}.txt \
    | ./scripts/deseq-output2gff ${ACCESSION} ${CHANGE_CUTOFF} \
				 > ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff
set +e # if egrep matches nothing
cat ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff \
    | egrep '; colour [23];' > ${DESEQ2}/changed${TAG}_${CHANGE_CUTOFF}.gff
set -e

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'

#! /bin/bash -x

# https://davetang.org/wiki/tiki-index.php?page=SAM#Clipped_alignment

samtools view -h data/03_star/aligned_0.bam \
    | ./scripts/sam2gff -p \
    | wc -l
#! /bin/bash

. doit-preamble.bash

DESEQ2=data/05_deseq2

# ------------------------------------------------------------------------
# Package the results
# ------------------------------------------------------------------------

RESULTS=$(date +%Y-%m-%d)-results
echo 1>&2 '# Creating '$RESULTS'...'

rm -rf ${RESULTS}
mkdir -p ${RESULTS}

#cp ${PROFILES}/*.profile ${RESULTS}
cp ${DESEQ2}/results_*.txt  ${RESULTS}

for f in ${DESEQ2}/*.gff ; do
    name=$(basename $f .gff)
    cat $f | ./scripts/split-gff -n -d ${RESULTS} $name
done

#cp ${STATS}/stats.txt ${RESULTS}/stats.txt
# ./scripts/tsv2xlsx -o ${RESULTS}/stats.xlsx ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
