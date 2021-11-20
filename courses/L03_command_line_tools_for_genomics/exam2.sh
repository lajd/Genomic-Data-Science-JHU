# Note: this file assumes that the gencommand_proj2_data.tar.gz file has been extracted into this directory

# Problem 1
## Number of alignments
NUM_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | wc -l)
echo "Number of alignments is $NUM_ALIGNMENTS"


# Problem 2
## Number unmapped alignments
NUM_UNMAPPED_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f7 | grep "*" | wc -l)
echo "Number of unmapped alignments is $NUM_UNMAPPED_ALIGNMENTS"

# Problem 3
## Number of alignments with a deletion
## Extract CIGAR string and check for deletions (D)
NUM_ALIGNMENTS_WITH_DELETION=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f6 | grep "D" | wc -l)
echo "Number of alignments with a deletion is $NUM_ALIGNMENTS_WITH_DELETION"

# Problem 4
## Number of alignemnts where read mate is mapped to the same chromosome
NUM_ALIGNMENTS_MATE_SAME_CHROMOSOME=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f7 | grep "=" | wc -l)
echo "Number of alignments with mate on same chromosome $NUM_ALIGNMENTS_MATE_SAME_CHROMOSOME"

# Problem 5
## Number of spliced alignments
NUM_SPLICED_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f6 | grep "N" | wc -l)
echo "Number of spliced alignments $NUM_SPLICED_ALIGNMENTS"


##
# Extract only the alignments in the range “Chr3:11,777,000-11,794,000”, corresponding to a
# locus of interest
##


# Sort and index to improve lookup
samtools sort gencommand_proj2_data/athal_wu_0_A.bam gencommand_proj2_data/athal_wu_0_A.sorted
samtools index gencommand_proj2_data/athal_wu_0_A.sorted.bam

# Problem 6
## How many alignments (Chr3:11,777,000-11,794,000) does the set contain
NUM_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | wc -l)
echo "Number of alignments in set $NUM_ALIGNMENTS"


# Problem 7
## How many alignments show that the read's mate is unmapped
NUM_UNMAPPED_MATES=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f7 | grep '*' | wc -l)
echo "Number of unmapped alignment mates in set $NUM_UNMAPPED_MATES"

# Problem 8
## How many alignments contain a deletion
NUM_ALIGNMENTS_WITH_DELETION=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f6 | grep 'D' | wc -l)
echo "Number of alignment with deletion in set $NUM_ALIGNMENTS_WITH_DELETION"

# Problem 9
## How many alignments show read mate mapped to the same chromosome
NUM_ALIGNMENTS_MATE_SAME_CHROMOSOME=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f7 | grep '=' | wc -l)
echo "Number of alignment with mate on the same chromosome $NUM_ALIGNMENTS_MATE_SAME_CHROMOSOME"

# Problem 10
## How many alignments are spliced
NUM_SPLICED_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f6 | grep 'N' | wc -l)
echo "Number of spliced alignments in set $NUM_SPLICED_ALIGNMENTS"


# Problem 11
## How many sequences are in the BAM file

echo "Number of sequences in BAM file is $( samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | grep -i "@sq" | wc -l)"

# Problem 12
## Length of first sequence in BAM file
echo "Length of the first sequence in BAM file is $( samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | grep -i "@sq" | head -n 1 | cut -f3 | grep -o -E '[0-9]+')"

# Problem 13
## Name of program used
echo "Name of program used is $( samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | grep -i "@PG"  | head -n 1 | cut -f2)"

# Problem 14
## What is the read identifier (name) for the first alignment?
echo "Read ID of first alignmentis $( samtools view gencommand_proj2_data/athal_wu_0_A.bam | head -n 1 | cut -f1 )"

# Problem 15
## What is the start position of this read’s mate on the genome? Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped.
samtools view gencommand_proj2_data/athal_wu_0_A.bam | head -n 1
# Mate is on same chromosome (chr3), and starts at pos 11700332
# chr3:11700332

# Problem 16
## How many overlaps (each overlap is reported on one line) are reported?
### get intersection of exons (gtf format) with alignments
echo "Num overlapping exon/alignments $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | wc -l)"

# Problem 17
## Base lengths (contained in col 16) greater than 10
echo "Num overlapping exon/alignments greater than 10 bases long is: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | cut -f16 | awk '$1 > 9' | wc -l)"

# Problem 18
## how many alignments overlap the annotation
echo "How many alignments overlap the annotations: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | wc -l)"


# Problem 19
## How many exons have reads mapped to them
## Identify exons as being the same if they have the same start location (f4)
echo "How many (unique) exons have reads mapped to them: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam  |  cut -f4 | sort -u | wc -l)"


# Problem 20
## If you were to convert the transcript annotations in the file “athal_wu_0_A_annot.gtf” into BED format, how many BED records would be generated?
## There would be a BED record for each unique gene
echo "Number of unique genes with alignments is: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | cut -f9 | cut -d " " -f4 | sort -u | wc -l)"
