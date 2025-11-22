#!/bin/bash -l
#$ -P vtrs
#$ -N BRB.Seq.stimulated
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -l h_rt=24:00:00


module load star/2.7.9a

STAR --runMode genomeGenerate \
	--genomeDir /rprojectnb/cancergrp/brb/STAR_Index \
	--genomeFastaFiles /rprojectnb/cancergrp/brb/genomes/GRCh38.primary_assembly.genome.fa \
	/rprojectnb/cancergrp/brb/genomes/all_vTRs.fasta \
	--sjdbGTFfile /rprojectnb/cancergrp/brb/annotation_files/gencode.v46_vTRs.gtf \
	--sjdbOverhang 79
