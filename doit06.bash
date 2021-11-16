#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 6. Make count tables
# ------------------------------------------------------------------------

echo 1>&2 '# Making count tables'

rm -rf ${COUNTS}
mkdir -p ${COUNTS}

init_FEATURECOUNTS_ARGS

FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level

fgrep $'\t'gene$'\t' ${INPUTS}/annotation-host.gtf \
      | ./scripts/sanitize-gtf-for-featureCounts \
	    > ${COUNTS}/annotation-host.gtf
fgrep $'\t'gene$'\t' ${INPUTS}/annotation-bacteria.gtf \
      | ./scripts/sanitize-gtf-for-featureCounts \
	    > ${COUNTS}/annotation-bacteria.gtf
cat ${COUNTS}/annotation-host.gtf ${COUNTS}/annotation-bacteria.gtf \
    > ${COUNTS}/annotation-both.gtf


(
    for i in $SAMPLES_INDICES ; do

	echo featureCounts $FEATURECOUNTS_ARGS \
    	     -a ${COUNTS}/annotation-both.gtf \
    	     -o ${COUNTS}/both_$i.txt \
    	     ${STAR}/aligned_$i.bam

	echo featureCounts $FEATURECOUNTS_ARGS \
	     -a ${COUNTS}/annotation-host.gtf \
	     -o ${COUNTS}/host_$i.txt \
	     ${STAR}/aligned_$i.bam

	echo featureCounts $FEATURECOUNTS_ARGS \
	     -a ${COUNTS}/annotation-bacteria.gtf \
	     -o ${COUNTS}/bacteria_$i.txt \
	     ${STAR}/aligned_$i.bam

    done
) | run_commands

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
