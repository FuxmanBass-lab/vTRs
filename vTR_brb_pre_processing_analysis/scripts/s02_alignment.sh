#!/bin/bash

# List of all samples
samples=(
    "Batch1.stimulated"
    "Batch1.unstimulated"
    "Batch2.stimulated"
    "Batch2.unstimulated"
    "Batch3.stimulated"
    "Batch3.unstimulated"
    "Batch4.stimulated"
    "Batch4.unstimulated"
    "Batch5.stimulated"
    "Batch5.unstimulated"
    "Batch6.stimulated"
    "Batch6.unstimulated"
)

# Iterate over each sample and create a separate job script
for sample in "${samples[@]}"; do
    # Create a temporary script file for each sample
    cat << EOF > temp_script_${sample}.sh
#!/bin/bash -l
#$ -P vtrs
#$ -N step06_brb_alignment_${sample}
#$ -l mem_per_core=16G
#$ -pe omp 8
#$ -l h_rt=12:00:00

module load star

STAR --soloType CB_UMI_Simple \
    --outSAMtype BAM SortedByCoordinate \
    --clipAdapterType CellRanger4 \
    --outSAMunmapped Within \
    --limitOutSJcollapsed 2000000 \
    --outReadsUnmapped Fastx \
    --soloFeatures Gene \
    --outSAMattributes NH HI AS nM CR CY UR UY CB UB GX GN sS sQ sM \
    --soloStrand Forward \
    --runThreadN 4 \
    --sjdbGTFfile /rprojectnb/cancergrp/brb/annotation_files/gencode.v46_vTRs.gtf \
    --genomeDir /rprojectnb/cancergrp/brb/STAR_Index \
    --soloCBstart 1 \
    --soloCBlen 14 \
    --soloUMIstart 15 \
    --soloUMIlen 14 \
    --soloCellFilter None \
    --soloCBwhitelist /rprojectnb/cancergrp/brb/raw_data/barcodes.txt \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_Directional \
    --readFilesCommand zcat \
    --outFileNamePrefix /rprojectnb/cancergrp/brb/intermediate_outputs/${sample}/out. \
    --readFilesIn /projectnb/vtrs/BRBSEQ_fastq/${sample}.2.fastq.gz /projectnb/vtrs/BRBSEQ_fastq/${sample}.1.fastq.gz
EOF

    # Submit the job to the queue
    qsub temp_script_${sample}.sh

done
