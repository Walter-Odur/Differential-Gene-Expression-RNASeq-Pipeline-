#!/bin/bash

# RNASeq Analysis

SECONDS=0

ls ./reads/*fastq |sed 's/_[^_]*.fastq$//g' |sort -u |xargs -i{} basename {} > sample.txt

# make necessary directories
mkdir -p fastqc trim starIdx bams counts


# Obtain Sample IDs
ids=($(cat sample.txt))

#------
# STEP1: Quality Control
#------
for id in "${ids[@]}"; do
        fastqc \
                reads/${id}*.fastq \
                -o ./fastqc
done

#------
# STEP2: Trimming
#------
for id in "${ids[@]}"; do
        trimmomatic PE \
                -threads 4 \
                -phred33 \
                ./reads/${id}_1.fastq ./reads/${id}_2.fastq \
                ./trim/${id}_1_trimmed.fastq ./trim/${id}_1_unpaired.fastq \
                ./trim/${id}_2_trimmed.fastq ./trim/${id}_2_unpaired.fastq \
                TRAILING:10
done

#------
# STEP3: Quality Control
#------
for id in "${ids[@]}"; do
        fastqc \
                ./trim/${id}*_trimmed.fastq \
                -o ./trim
done


: << 'COMM'
#------
# STEP4: Index The Reference    ----- If STAR indexes are unavailable
#------
STAR \
        --runThreadN 4 \
        --runMode genomeGenerate \
        --genomeDir ./starIdx \
        --genomeFastaFiles ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --sjdbGTFfile ./ref/Homo_sapiens.GRCh38.112.gtf \
        --sjdbOverhang 150
COMM



#------
# STEP5: Alignment
#------
for id in ${ids[@]}; do
        STAR \
                --runThreadN 4 \
                --genomeDir ./starIdx \
                --readFilesIn ./trim/${id}_1_trimmed.fastq ./trim/${id}_2_trimmed.fastq \
                --outFileNamePrefix ./bams/${id}_ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard
done

#------
# STEP6: Index Sorted bam files
#------
for id in ${ids[@]}; do
        samtools index \
                bams/${id}_*.out.bam
done

#------
# STEP7: Quantification
#------
bams=(bams/*.out.bam)
annotFile="ref/Homo_sapiens.GRCh38.112.gtf"
output="counts/feature_counts.txt"

featureCounts -T 4 -a ${annotFile} -o ${output} "${bams[@]}"


# Calculate duration

duration=${SECONDS}

echo "$((${duration} / 3600)) hrs: $((($duration % 3600) / 60)) mins: $(($duration % 60)) secs"