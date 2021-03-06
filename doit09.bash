#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Package the results
# ------------------------------------------------------------------------

echo 1>&2 fixme: need to check for missing edits
exit 1

RESULTS=$(date +%Y-%m-%d)-results
echo 1>&2 '# Creating '$RESULTS'...'

rm -rf ${RESULTS}
mkdir -p ${RESULTS}

cp ${PROFILES}/*.profile ${RESULTS}
cp ${DESEQ2}/results_*.txt  ${RESULTS}

for f in ${DESEQ2}/*.gff ; do
    name=$(basename $f .gff)
    cat $f | ./scripts/split-gff -n -d ${RESULTS} $name
done

 cp ${STATS}/stats.txt ${RESULTS}/stats.txt
# ./scripts/tsv2xlsx -o ${RESULTS}/stats.xlsx ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
