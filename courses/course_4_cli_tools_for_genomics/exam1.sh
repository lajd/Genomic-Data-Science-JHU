# Note: this file assumes that the gencommand_proj1_data.tar.gz file has been extracted into this directory

# Problem 1
## Number of chromosomes in the genome
NUM_CHROMOSOMES=$(grep ">" gencommand_proj1_data/apple.genome | wc -l)
echo "Number of chromosomes is $NUM_CHROMOSOMES"

# Problem 2
## Number of genes
NUM_GENES=$( cut -f1 gencommand_proj1_data/apple.genes | uniq | wc -l)
echo "Number of genes is $NUM_GENES"

# Problem 3
## Number of transcript variants
NUM_TRANSCRIPT_VARIANTS=$( cut -f2 gencommand_proj1_data/apple.genes | uniq | wc -l)
echo "Number of transcript variants is $NUM_TRANSCRIPT_VARIANTS"

# Problem 4
## Number of single splice variant genes
NUM_SSV_GENES=$( cut -f1 gencommand_proj1_data/apple.genes | uniq -c | grep "1 " | wc -l)
echo "Number of single splice variant genes is $NUM_SSV_GENES"

# Problem 5
## Number of double splice variant genes
NUM_DSV_GENES=$( cut -f1 gencommand_proj1_data/apple.genes | uniq -c | grep "2 " | wc -l)
echo "Number of double splice variant genes is $NUM_DSV_GENES"


# Problem 6
## Number of genes on + strand
NUM_PLUS_GENES=$(cut -f4 -f1 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "+" | wc -l)
echo "Number of genes on + strand is $NUM_PLUS_GENES"

# Problem 7
## Number of genes on - strand
NUM_NEG_GENES=$(cut -f4 -f1 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "-" | wc -l)
echo "Number of genes on - strand is $NUM_NEG_GENES"

# Problem 8
## Number of genes on chr1
NUM_GENES_CHR1=$(cut -f3 -f1 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "chr1" | wc -l)
echo "Number of genes on chr1 $NUM_GENES_CHR1"

# Problem 9
## Number of genes on chr2
NUM_GENES_CHR2=$(cut -f3 -f1 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "chr2" | wc -l)
echo "Number of genes on chr2 $NUM_GENES_CHR2"


# Problem 10
## Number of genes on chr3
NUM_GENES_CHR3=$(cut -f3 -f1 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "chr3" | wc -l)
echo "Number of genes on chr3 $NUM_GENES_CHR3"


# Problem 11
## Number of transcripts on chr1
NUM_TRANS_CHR1=$(cut -f3 -f2 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "chr1" | wc -l)
echo "Number of transcripts on chr1 $NUM_TRANS_CHR1"


# Problem 12
## Number of transcripts on chr2
NUM_TRANS_CHR2=$(cut -f3 -f2 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "chr2" | wc -l)
echo "Number of transcripts on chr2 $NUM_TRANS_CHR2"

# Problem 13
## Number of transcripts on chr3
NUM_TRANS_CHR3=$(cut -f3 -f2 gencommand_proj1_data/apple.genes | sort -u -k1,1 | grep "chr3" | wc -l)
echo "Number of transcripts on chr3 $NUM_TRANS_CHR3"


# Problem 14
## Number of genes in common between A and B
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
NUM_COMMON_GENES_A_B=$(comm -1 -2 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique | wc -l)
echo "Number of common genes between A and B is $NUM_COMMON_GENES_A_B"

# Problem 15
## Number of genes specific to A
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
NUM_SPECIFIC_GENES_A=$(comm -2 -3 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique | wc -l)
echo "Number of specific genes to A is $NUM_SPECIFIC_GENES_A"



# Problem 16
## Number of genes specific to B
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
NUM_SPECIFIC_GENES_B=$(comm -1 -3 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique | wc -l)
echo "Number of specific genes to A is $NUM_SPECIFIC_GENES_B"


# Problem 16
## Number of genes common to A, B, C
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
comm -1 -2 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique > /tmp/apple.condition_common_to_A_B.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionC | sort -u > /tmp/apple.conditionC.sorted.unique

NUM_COMMON_GENES_A_B_C=$(comm -1 -2 /tmp/apple.condition_common_to_A_B.sorted.unique /tmp/apple.conditionC.sorted.unique | wc -l)
echo "Number of common genes to A, B and C are $NUM_COMMON_GENES_A_B_C"


