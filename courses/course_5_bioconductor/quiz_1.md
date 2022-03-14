# Quiz 1

## Setup
```r
library(AnnotationHub)
ah = AnnotationHub()
```

## Question 1

```r
cpgislands = query(subset(ah, species == "Homo sapiens"), "CpG Islands")

# Take most recent version (hg19)
cpgislands$genome

hg19 = cpgislands[[1]]

# Keep only standard chromosomes
standard_chromosomes = keepStandardChromosomes(hg19, pruning.mode = "coarse")

# Remove sex chromosomes
cpg_autosomes = dropSeqlevels(standard_chromosomes, c("chrX", "chrY", "chrM"), pruning.mode = "coarse")

num_islands_autosomes = length(cpg_autosomes)
# num_islands_autosomes=26641
```


## Question 2
```r
# Keep only standard chromosomes
standard_chromosomes = keepStandardChromosomes(hg19, pruning.mode = "coarse")
chr4 = keepSeqlevels(standard_chromosomes, "chr4", pruning.mode = "coarse")
num_islands_chr4 = length(chr4)
# num_islands_chr4=1031
```

## Question 3
```r

## Further filter for narrow peaks 
h1_H3k4me3_granges = subset(ah, species == "Homo sapiens") %>%
    query(., c("H3k4me3", "H1", "Granges", "Epigenome", "Narrow"))

# Arbitrarily use the first one
h1_H3k4me3_granges_narrow_peak = h1_H3k4me3_granges[[1]]

# Only take the autosomes
H3k4me3_autosomes = keepStandardChromosomes(h1_H3k4me3_granges_narrow_peak, pruning.mode = "coarse") %>%
    dropSeqlevels(., c("chrX", "chrY"), pruning.mode = "coarse")

# Sum the coverage of the widths
num_bps_covered = sum(sum(coverage(reduce(H3k4me3_autosomes))))
#num_bps_covered=41135164
```

## Question 4

```r
# Further filter for narrow peaks 
h1_H3K27me3_granges = subset(ah, species == "Homo sapiens") %>%
    query(., c("H3K27me3", "H1", "Granges", "Epigenome", "Narrow"))

# Arbitrarily use the first one
h1_H3K27me3_granges_narrow_peak = h1_H3K27me3_granges[[1]]

# Only take the autosomes
H3K27me3_autosomes = keepStandardChromosomes(h1_H3K27me3_granges_narrow_peak, pruning.mode = "coarse") %>%
    dropSeqlevels(., c("chrX", "chrY"), pruning.mode = "coarse")

# Get mean
mean_signal_value = mean(H3K27me3_autosomes$signalValue)
# mean_signal_value=4.770728
```


## Question 5
```r
# Find bivalent region by taking the intersection of the H3K4me3 and H3K27me3 regions
bivalent_intervals = intersect(H3K27me3_autosomes, H3k4me3_autosomes)
num_bivalent_bps = sum(width(bivalent_intervals))
# num_bivalent_bps=10289096
```

## Question 6
```r
# Find the overlapping regions between the bivalent base pairs and the cpg autosomes
bivalent_cpg_overlaps = findOverlaps(bivalent_intervals, cpg_autosomes)

unique_overlap_hits = unique(queryHits(bivalent_cpg_overlaps))

fraction_bivalent_cpg_overlaps = length(unique_overlap_hits)/length(bivalent_intervals)
# fraction_bivalent_cpg_overlaps = 0.5383644
```

## Question 7
```r
# Find the overlapping bases between the bivalent base pairs and the cpg autosomes
bivalent_cpg_intersect = intersect(bivalent_intervals, cpg_autosomes)

fraction_bivalent_cpg_intersect = sum(width(reduce(bivalent_cpg_intersect))) / sum(width(cpg_autosomes))

# fraction_bivalent_cpg_intersect = 0.241688
```

## Question 8
```r

# Resize the cpg islands to add 10kb to each end of the region
cpg_autosomes_10kb = resize(cpg_autosomes, fix = "center", width = width(cpg_autosomes) + 20000)
# Intersect the new regions with the bivalent intervals
bivalent_cpg_10kb_intersect = intersect(bivalent_intervals, cpg_autosomes_10kb)

intersecting_bps = sum(width(bivalent_cpg_10kb_intersect))
# intersecting_bps = 9782086
```

## Question 9
```r
# calculate human genome size
num_genome_bps = sum(as.numeric(seqlengths(cpg_autosomes)))
# Calculate the number of bps in cpg intervals
cpg_interval_size = sum(width(cpg_autosomes))
frac_genome_in_cpg_island = cpg_interval_size / num_genome_bps
# frac_genome_in_cpg_island =  0.007047481
```


## Question 10
```r
# calculate human genome size
odds_matrix = matrix(0, ncol=2, nrow=2)

colnames(odds_matrix) = c("in", "out")
rownames(odds_matrix) = c("in", "out")


odds_matrix[1,1] = sum(width(bivalent_cpg_intersect))
odds_matrix[2,1] = sum(width(setdiff(cpg_autosomes, bivalent_intervals)))
odds_matrix[1,2] = sum(width(setdiff(bivalent_intervals, cpg_autosomes)))
odds_matrix[2,2] = num_genome_bps - sum(odds_matrix)

odds_ratio = odds_matrix[1,1] * odds_matrix[2,2] / (odds_matrix[2,1] * odds_matrix[1,2])
# odds_ratio=169.0962
```
