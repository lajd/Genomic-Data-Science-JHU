# Note: this file assumes that the gencommand_proj3_data.tar.gz file has been extracted into this directory

# Problem 1
## Number of sequences in the genome
NUM_SEQUENCES=$(cat gencommand_proj3_data/wu_0.v7.fas | grep ">" | wc -l)
echo "Number of sequences is $NUM_SEQUENCES"


# Problem 2
## Name of third sequence
THIRD_SEQUENCE=$(cat gencommand_proj3_data/wu_0.v7.fas | grep ">" | head -n 3 | tail -n 1)
echo "Name of third sequence is $THIRD_SEQUENCE"

# Problem 3
## Name of last sequence
LAST_SEQUENCE=$(cat gencommand_proj3_data/wu_0.v7.fas | grep ">" | tail -n 1)
echo "Name of last sequence is $LAST_SEQUENCE"

# Problem 4
## Generate the Bowtie2 index
mkdir -p gencommand_proj3_data/wu_0_bt2/
bowtie2-build gencommand_proj3_data/wu_0.v7.fas gencommand_proj3_data/wu_0_bt2/wu_0
## Number of bt2 files is
NUM_BT2_FILES=$(ls gencommand_proj3_data/wu_0_bt2/wu_0*.bt2 | wc -l)
echo "Name of bt2 index files is $NUM_BT2_FILES"

# Problem 5
## Extension is "bt2"


# Problem 6
## Number of reads in original fastq file
NUM_LINES=$(cat gencommand_proj3_data/wu_0_A_wgs.fastq | wc -l)
NUM_FASTQ_READS=$((NUM_LINES/4))
echo "Name of fastq reads is $NUM_FASTQ_READS"


# Problem 7/9/11
## Number of matched & mapped alignments
## Align fastq reads to the created bt2 index, and output as gencommand_proj3_data/wu_0.bt2.sam
bowtie2 -p 4 -x gencommand_proj3_data/wu_0_bt2/wu_0 gencommand_proj3_data/wu_0_A_wgs.fastq -S gencommand_proj3_data/wu_0.bt2.sam
# Output:
#147354 reads; of these:
#  147354 (100.00%) were unpaired; of these:
#    9635 (6.54%) aligned 0 times
#    93780 (63.64%) aligned exactly 1 time
#    43939 (29.82%) aligned >1 times
#93.46% overall alignment rate

echo "Number of matched & mapped alignments is 147354-9635=137719"
echo "Number of matched alignments is 147354-9635=137719"
echo "Number of multi-alignemnts 43939"


# Problem 8/10/12
## Number of matched & mapped alignments in local match setting
## Align fastq reads to the created bt2 index, using local alignment,
## and output as gencommand_proj3_data/wu_0.bt2.local.sam
bowtie2 -p 4 -x gencommand_proj3_data/wu_0_bt2/wu_0 gencommand_proj3_data/wu_0_A_wgs.fastq -S gencommand_proj3_data/wu_0.bt2.local.sam  --local
# Output:
#147354 reads; of these:
#  147354 (100.00%) were unpaired; of these:
#    6310 (4.28%) aligned 0 times
#    84939 (57.64%) aligned exactly 1 time
#    56105 (38.07%) aligned >1 times
#95.72% overall alignment rate
echo "Number of matched & mapped (local) alignments is 147354-6310=141044"
echo "Number of matched (local) alignments is 147354-6310=141044"
echo "Number of multi-alignemnts 56105"


# Problem 13
## How many alignments contained insertions and/or deletions, in the scenario in Question 7
NUM_ALIGNMENTS_WITH_INDELS=$(cat gencommand_proj3_data/wu_0.bt2.sam | cut -f6 | grep "I\|D" | wc -l)
echo "Number of alignments with indels is: $NUM_ALIGNMENTS_WITH_INDELS"

# Problem 14
## How many alignments contained insertions and/or deletions, in the scenario in Question 7
NUM_ALIGNMENTS_WITH_INDELS=$(cat gencommand_proj3_data/wu_0.bt2.local.sam | cut -f6 | grep "I\|D" | wc -l)
echo "Number of alignments with indels is: $NUM_ALIGNMENTS_WITH_INDELS"


######
#Compile candidate sites of variation using SAMtools mpileup for further evaluation with BCFtools.
# Provide the reference fasta genome and use the option “-uv” to generate the output in uncompressed
# VCF format for easy examination.
#
# Align and output as bam
samtools view -bT gencommand_proj3_data/wu_0.v7.fas gencommand_proj3_data/wu_0.bt2.sam > gencommand_proj3_data/wu_0.bt2.bam
# Sort the bam
samtools sort gencommand_proj3_data/wu_0.bt2.bam gencommand_proj3_data/wu_0.bt2.sorted.bam
# Index the bam
samtools index gencommand_proj3_data/wu_0.bt2.sorted.bam
# Create VCF format using mpileup of the reference against the alignments
samtools mpileup -uv -f gencommand_proj3_data/wu_0.v7.fas gencommand_proj3_data/wu_0.bt2.sorted.bam > gencommand_proj3_data/wu_0.vcf


# Problem 15
# Count the variations which occur on Chr3
NUM_VARIATIONS_ON_CHR3=$(cat gencommand_proj3_data/wu_0.vcf | cut -f1 | grep "Chr3" | wc -l)
echo "Number of alignments on Chr3 is: $NUM_VARIATIONS_ON_CHR3"

# Problem 16
#How many entries have ‘A’ as the corresponding genome letter?
NUM_VARIATIONS_WITH_REF_A=$(cat gencommand_proj3_data/wu_0.vcf | cut -f4 | grep "A" | wc -l)
echo "Number of variations with reference=A is: $NUM_VARIATIONS_WITH_REF_A"

# Problem 17
NUM_ENTRIES_WITH_READ_DEPTH_20=$(cat gencommand_proj3_data/wu_0.vcf | cut -f8 | grep "DP=20;" | wc -l)
echo "Number of entries with read depth 20 is: $NUM_ENTRIES_WITH_READ_DEPTH_20"

# Problem 18
NUM_ENTRIES_WITH_INDELS=$(cat gencommand_proj3_data/wu_0.vcf | cut -f8 | grep "INDEL;" | wc -l)
echo "Number of entries with indels is: $NUM_ENTRIES_WITH_INDELS"


# Problem 19
# Number of entries on chr2 with position 175672
NUM_EXTRIES_ON_CHR1_WITH_POS=$(cat gencommand_proj3_data/wu_0.vcf | cut -f1 -f2 | grep "Chr1\t175672$" | wc -l)
echo "Number of entries on chr1 with position 175672 is: $NUM_EXTRIES_ON_CHR1_WITH_POS"


# Problem 20
## How many variants are called on Chr3

# First create BCF format
samtools mpileup -g -f gencommand_proj3_data/wu_0.v7.fas gencommand_proj3_data/wu_0.bt2.sorted.bam > gencommand_proj3_data/wu_0.bcf
## Call BCFtools on the created VCF file (variants only)
bcftools call -m -v -O z -o gencommand_proj3_data/wu_0.vcf.gz gencommand_proj3_data/wu_0.bcf
# Use zcat to view
NUM_VARIANTS_ON_CHR3=$(gzcat gencommand_proj3_data/wu_0.vcf.gz | cut -f1 | grep "Chr3$" | wc -l)
echo "Number of entries on chr1 with position 175672 is: $NUM_VARIANTS_ON_CHR3"


# Problem 21
## A->T SNP
# Use zcat to view
NUM_VARIANTS_A_TO_T=$(gzcat gencommand_proj3_data/wu_0.vcf.gz | cut -f4 -f5 | grep "^A\tT$" | wc -l)
echo "Number of variants A to T is: $NUM_VARIANTS_A_TO_T"

# Problem 22
## How many entries are indels
# Use zcat to view
NUM_ENTRIES_WITH_READ_DEPTH_20=$(gzcat gencommand_proj3_data/wu_0.vcf.gz | cut -f8 | grep "DP=20" | wc -l)
echo "Number entries with read depth 20 is: $NUM_ENTRIES_WITH_READ_DEPTH_20"


# Problem 23
## What type of variant (i.e., SNP or INDEL) is called at position 11937923 on Chr3?
gzcat gencommand_proj3_data/wu_0.vcf.gz  | grep "Chr3" | grep "11937923"
echo "Type of variant is: SNP"
