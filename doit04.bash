#! /bin/bash -x

# https://davetang.org/wiki/tiki-index.php?page=SAM#Clipped_alignment

samtools view -h data/03_star/aligned_0.bam \
    | ./scripts/sam2gff -p \
    | wc -l
