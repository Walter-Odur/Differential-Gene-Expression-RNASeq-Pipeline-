#!/bin/bash


echo "
Steps to generate featureCount files:
WARNING!
PLEASE carefully go through the README file before executing this script!
"

SECONDS=0

# create necessary directories
mkdir -p fastqc trim starIdx bams counts

# create sample ID file
ls ./reads/* fastq |sed 's/_[ˆ_]*.fastq$//g' |\
    sort -u |xargs -i{} basename {} > samples.txt

# assign a variable name to sample IDs
ids=$(cat samples.txt)


#----
# STEP1: QUALITY CONTROL
#----

# run FASTQC
echo "running fastqc..."
for id in ${ids[@]}
do
        fastqc \
                reads/${id}* .fastq \
                -o fastqc
done


# trim poor quqlity reads
echo "trimming poor quality reads..."
for id in ${ids[@]}
do
        trimmomatic \
                PE \
                -threads 4 \
                -phred33 \
                reads/${id}_1.fastq reads/${id}_2.fastq \
                trim/${id}_1_trimmed.fastq trim/${id}_1_unpaired.fastq \
                trim/${id}_2_trimmed.fastq trim/${id}_2_unpaired.fastq \
                TRAILING:10
done


#----
# STEP2: INDEX REFERENCE GENOME
#----
if ! ls starIdx 1>/dev/null 2>&1
then
        echo "Indexing the reference..."
        STAR \
                --runThreadN 4 \
                --runMode genomeGenerate \
                --genomeDir starIdx \
                --genomeFastaFiles ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                --sjdbGTFfile ref/Homo_sapiens.GRCh38.113.gtf \
                --sjdbOverhang 50 # readLength - 1 [Length for splice junction database]
else
        echo "reference indexes found!"
fi


#----
# STEP3: READ ALIGNMENT
#----
echo "read alignment..."
for id in ${ids[@]}
do
        STAR \
                --runThreadN 4 \
                --genomeDir starIdx \
                --readFilesIn trim/${id}_1_trimmed.fastq trim/${id}_2_trimmed.fastq \
                --outFileNamePrefix bams/${id}_ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard
done


# Index BAM files
echo "indexing BAM files..."
for id in ${ids[@]}
do
        samtools index \
                bams/${id}_*.out.bam
done


#----
# STEP4: FEATURE QUANTIFICATION
#----
echo "feature quantification..."
bams=(bams/*.out.bam)
annotFile="ref/Homo_sapiens.GRCh38.112.gtf"
output="counts/feature_counts.txt"

featureCounts \
                -T 4 \
                -a ${annotFile} \
                -o ${output} "${bams[@]}"


## Report the process duration
duration=${SECONDS}
echo "
The whole process took:
$(($duration / 3600)) hours: $((($duration % 3600) / 60)) minutes: $(($duration % 60)) seconds
"
