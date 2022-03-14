# Quiz 2

# Setup

```r
library(BSgenome)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("BSgenome")

library(AnnotationHub)
ah=AnnotationHub()

available.genomes()

install("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")

seqlengths(BSgenome.Hsapiens.UCSC.hg19)
seqnames(BSgenome.Hsapiens.UCSC.hg19)

bsgenome_hg19_chr22 = BSgenome.Hsapiens.UCSC.hg19$chr22

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

transcripts_db = transcripts(txdb)
genes_db = genes(txdb)
transcript_lengths_db = transcriptLengths(txdb, with.cds_len = TRUE)

# Use chr22 length = 51304566 (wikipedia)
chr22_length = 51304566
chr22_gr = GRanges(seqnames="chr22", ranges=IRanges(start=1, end=chr22_length))
```

# Q1

```r
ah_subset = query(ah,c("H3K27me3","narrowPeak", "E003"))
ah_29892 = ah_subset[["AH29892"]]

atcg_sum = letterFrequency(bsgenome_hg19_chr22, "ATCG", as.prob=FALSE)
gc_sum = letterFrequency(bsgenome_hg19_chr22, "GC", as.prob=FALSE)
gc_frac = gc_sum / atcg_sum
gc_frac
# gc_frac = 0.4798807
```

# Q2
```r
ah_subset = query(ah,c("H3K27me3","narrowPeak", "E003"))
ah_29892 = ah_subset[["AH29892"]]
ah_29892_chr22 = keepSeqlevels(ah_29892, c("chr22"), pruning.mode="tidy")

# Trim
trim(ah_29892_chr22)

ah_29892_chr22_view = Views(Hsapiens, ah_29892_chr22)

letterFrequency(ah_29892_chr22_view,"GC", as.prob=TRUE)

# Avg GC content on peaks
gc_content_for_peak = mean(letterFrequency(ah_29892_chr22_view, "GC", as.prob=TRUE))
gc_content_for_peak
# gc_content_for_peak = 0.528866
```

# Q3
```r
ah_29892_chr22_signal_value = ah_29892_chr22$signalValue
# Create df aligning GC content with signal value

correlation_gc_signal = cor(letterFrequency(ah_29892_chr22_view, "GC", as.prob=TRUE), ah_29892_chr22_signal_value)

correlation_gc_signal
# correlation_gc_signal = 0.004467924
````

# Q4
```r
ah_chip_subset = query(ah, c("fc.signal", "H3K27me3", "E003"))
# Use AH32033 in ah_chip_subset
ah_32033_chip_bw = ah_chip_subset[["AH32033"]]

# import from ah_29892_chr22
ah_32033_chip = import(ah_32033_chip_bw, which=ah_29892_chr22)

# Only import the subset defined by chr22
ah32033_chip_chr22 = import(ah_32033_chip_bw, which=chr22_gr, as="Rle")$chr22

# Compute correlation
# Get fc.signal across the same narrow-peak regions defined by ah_29892_chr22
fc_signal_across_narrow_peak_view = Views(ah32033_chip_chr22, start=start(ah_29892_chr22), end=end(ah_29892_chr22))
correlation_fc_signal_np_signal = cor(mean(fc_signal_across_narrow_peak_view), ah_29892_chr22_signal_value)

correlation_fc_signal_np_signal
# correlation_fc_signal_np_signal = 0.9149614
```

# Q5
```r
num_bases_ge_1 = sum(ah32033_chip_chr22 >= 1)
num_bases_ge_1
# num_bases_ge_1=10914671
```

# Q6
```r
ah_e055_subset = query(ah,c("H3K27me3", "fc.signal", "E055"))
ah_32470_bw = ah_e055_subset[["AH32470"]]
ah_32470_chip_chr22 = import(ah_32470_bw, which=chr22_gr, as="Rle")$chr22


ah_32470_chip_chr22_ge_2 = (ah_32470_chip_chr22 >= 2)
ah32033_chip_chr22_le_05 = (ah32033_chip_chr22 <= 0.5)

ah_32470_chip_chr22_ge_2_gr = as(ah_32470_chip_chr22_ge_2, "IRanges")
ah32033_chip_chr22_le_05_gr = as(ah32033_chip_chr22_le_05, "IRanges")

ir_intersect = intersect(ah_32470_chip_chr22_ge_2_gr, ah32033_chip_chr22_le_05_gr)
bp_intersect = sum(width(ir_intersect))

bp_intersect
# bp_intersect = 1869937
```

# Q7

```r
# Query for CpG islands
CpG_islands_subset = query(ah, c("hg19", "CpG"))
CpG_islands = CpG_islands_subset[["AH5086"]]
CpG_islands_chr22 = keepSeqlevels(CpG_islands, c("chr22"), pruning.mode="tidy")

# Obtain the sequence data
cpg_island_view = Views(bsgenome_hg19_chr22, start=start(CpG_islands_chr22), end=end(CpG_islands_chr22))

# Expected frequency
CpG_G_count = letterFrequency(cpg_island_view,"G")
CpG_C_count = letterFrequency(cpg_island_view,"C")
expected_cpg_frequency = (CpG_G_count * CpG_C_count) / width(CpG_islands_chr22)

# Observed frequency
observed_dinucleotide_freq = dinucleotideFrequency(cpg_island_view)
observed_cpg_dinucleotide_freq = observed_dinucleotide_freq[,"CG"]

# Observed / expected
ratio_observed_expected_cpg = mean(observed_cpg_dinucleotide_freq / expected_cpg_frequency)
ratio_observed_expected_cpg
# ratio_observed_expected_cpg=0.8340929
```

# Q8

```
TATA_dna_str = DNAString("TATAAA")

TATA_match_chr22 = matchPattern(TATA_dna_str, bsgenome_hg19_chr22) 
TATA_rev_compl_match_chr22 = matchPattern(reverseComplement(TATA_dna_str), bsgenome_hg19_chr22) 
num_hits = length(TATA_match_chr22)+length(TATA_rev_compl_match_chr22)

num_hits
# num_hits=27263

```

# Q9

```
# Don't ignore the strand (same strand)
chr_22_transcripts = keepSeqlevels(transcripts_db, c("chr22"), pruning.mode="tidy")
chr_22_promoter_regions = promoters(chr_22_transcripts, upstream = 900, downstream = 100)

chr_22_genes = keepSeqlevels(genes_db, c("chr22"), pruning.mode="tidy")

# Overlap of promoter regions with genes (which promoter regions have a gene)
promoters_containing_genes = findOverlaps(chr_22_promoter_regions, chr_22_genes)

# For each of the regions containing a gene, identify whether the TATA box is contained
promoter_containing_gene_region_indices = unique(queryHits(promoters_containing_genes))

get_num_TATA_boxes <- function(idx) {
    promoter_region_seqdata = DNAStringSet(Views(Hsapiens, chr_22_promoter_regions[idx]))
    return(vcountPattern(TATA_dna_str, promoter_region_seqdata, with.indels = FALSE))
}

TATA_box_count_in_regions = lapply(promoter_containing_gene_region_indices, test_fun)

sum_TATA_boxes = Reduce("+", TATA_box_count_in_regions)
# sum_TATA_boxes = 266
# No matching answer -- choose closest answer 252

```

# Q10

```
# Get the promoter regions 900bp upstream and 100bp downstream on chr22
chr_22_transcripts = keepSeqlevels(transcripts_db, c("chr22"), pruning.mode="tidy")
chr_22_promoters = promoters(chr_22_transcripts, upstream = 900, downstream = 100)
# Get the transcript IDs of these promoters
chr_22_promoters_tx_ids = mcols(chr_22_promoters)$tx_id

# Get chr22 transcripts that have 1 or more associated CDS regions
tx_lengths_db_ge_1_cds_  = transcript_lengths_db[transcript_lengths_db$cds_len > 0,]
# Get transcript IDs containing CDS 
tx_ids_ge_1_cds  = tx_lengths_db_ge_1_cds_$tx_id

# Transaction ID intersection mask (promoter regions with at least 1 CDS)
tx_id_intersect_mask = chr_22_promoters_tx_ids %in% tx_ids_ge_1_cds

# Find the chr22 positions which are covered by 2 or more promoter regions containing a CDS 
chr22_bp_intervals_with_ge_2_promoters_of_cds = coverage(chr_22_promoters[tx_id_intersect_mask]) >= 2

bp_part_of_multiple_promoters_of_cds = sum(sum(chr22_bp_intervals_with_ge_2_promoters_of_cds))

bp_part_of_multiple_promoters_of_cds
# bp_part_of_multiple_promoters_of_cds = 306920
```
