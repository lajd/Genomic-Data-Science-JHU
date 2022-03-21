# Quiz 3

# Setup

```r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

install("ShortRead")
library(ShortRead)
install("BioStrings")
library(Biostrings)
install("yeastRNASeq")
library(yeastRNASeq)

install("leeBamViews")
library(leeBamViews)

install("Rsamtools")
library("Rsamtools")

install("oligo")
library(oligo)
instal("GEOquery")
library(GEOquery)

install("limma")
library(limma)

instal("minfi")
library(minfi)
library(minfiData)

install("AnnotationHub")
library(AnnotationHub)
ah = AnnotationHub()

install("DESeq2")
library(DESeq2)
install("zebrafishRNASeq")
library(zebrafishRNASeq)

```

# Q1

```r
fastq_file_path = system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
yeastrna_fastq = readFastq(fastq_file_path)
yeast_rnaseq_reads = sread(yeastrna_fastq)

yeastmat = as(yeast_rnaseq_reads, "matrix")
yeast_col_5 = table(yeastmat[, 5])

frac_of_reads_col_5 = yeast_col_5["A"] / sum(yeast_col_5)
frac_of_reads_col_5
# frac_of_reads_col_5 = 0.363841
```

# Q2
```r
qvals = quality(yeastrna_fastq)
yeastmat_q = as(qvals, "matrix")

avg_qval_col5 = mean(yeastmat_q[, 5])
avg_qval_col5
# avg_qval_col5 = 28.93
```

# Q3

```r
bamFilePath = system.file("bam", "isowt5_13e.bam", package="leeBamViews")
bamFile = BamFile(bamFilePath)
isowt5_reads = scanBam(bamFile)

yeast_chr_13_gr = GRanges(seqnames = "Scchr13", ranges = IRanges(start = c(800000), end = c(801000)))

# Read all reads from bamFile
reads = scanBam(bamFile, param=ScanBamParam(which = yeast_chr_13_gr, what = scanBamWhat()))

read_pos = reads[[1]]$pos
num_duplicated_reads_by_pos = sum(table(read_pos)[table(read_pos) > 1])
num_duplicated_reads_by_pos
# num_duplicated_reads_by_pos = 129
```

# Q4

```r
full_bam_paths = list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
gr_nt <- GRanges(seqnames="Scchr13", ranges=IRanges(start = c(807762), end = c(808068)))


get_mean_read_length <- function(file_path) {
  bamFile = BamFile(file_path)
  bamReads = scanBam(bamFile)
  reads = scanBam(bamFile, param=ScanBamParam(which = gr_nt, what = scanBamWhat()))
  return(length(reads[[1]][[1]]))
}

mean_read_length_across_samples = mean(unlist(lapply(full_bam_paths, get_mean_read_length)))
mean_read_length_across_samples
# mean_read_length_across_samples=90.25
```

# Q5

```r
# get data
library("GEOquery")
getGEOSuppFiles("GSE38792")
list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")

# read data; make sure to load all data (OSA/control)
gse38792_celfiles <- list.files("GSE38792/CEL", full = TRUE)
gse38792_celdata <- read.celfiles(gse38792_celfiles)

normalized_gse38792_celdata = rma(gse38792_celdata)
sample_8149273_row_idx = match("8149273", rownames(normalized_gse38792_celdata))
sample_8149273_expression = exprs(normalized_gse38792_celdata[sample_8149273_row_idx,])
mean_expr_level = mean(sample_8149273_expression)
mean_expr_level
# mean_expr_level = 7.039974
```

# Q6

```r
# Rename columns to control/OSA and get pData
# This is directly from tutorial
filename <- sampleNames(normalized_gse38792_celdata)
pData(normalized_gse38792_celdata)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(normalized_gse38792_celdata) <- sampleNames
pData(normalized_gse38792_celdata)$group <- ifelse(grepl("^OSA", sampleNames(normalized_gse38792_celdata)), "OSA", "Control")

# Compute design matrix
gse38792_design_matrix = model.matrix(~factor(normalized_gse38792_celdata$group))

fit = lmFit(normalized_gse38792_celdata, gse38792_design)
fit = eBayes(fit)

fit_toptable = topTable(fit)
log_fc = fit_toptable$logFC
abs_log_foldchange = abs(log_fc[1])
abs_log_foldchange
# abs_log_foldchange = 0.7126484
```

# Q7

```r
de_genes <- subset(fit_toptable, adj.P.Val < 0.05)
de_genes
# de_genes = 0
# There are no differentially expressed genes
```


# Q8

```r

# Load RGsetEx
rgSetExp = preprocessFunnorm(RGsetEx)
# Get the OPenSea CpG islands
openSeaRGSetExp = rgSetExp[c(getIslandStatus(rgSetExp) == "OpenSea")]
openSeaRGSetBeta = getBeta(openSeaRGSetExp)

pData(rgSetExp)[["Sample_Name"]]
# "GroupA_3" "GroupA_2" "GroupB_3" "GroupB_1" "GroupA_1" "GroupB_2"
# Group A is columns 1, 2, 5
# Group B is columns 3, 4, 6
# Take the first group (A) to be the control

openSeaBetaNormal = openSeaRGSetBeta[, c(1,2,5)]
openSeaBetaCancer = openSeaRGSetBeta[, c(3,4,6)]

# Compare the mean beta values
mean_beta_diff = mean(openSeaBetaNormal)-mean(openSeaBetaCancer)
mean_beta_diff
# mean_beta_diff=0.08863657
```


# Q9

```r
ah_caco2_subset = query(subset(ah, species == "Homo sapiens"), c("Caco2", "Dna", "Awg", "narrowPeak", "ENCODE"))
ah_caco2 = ah_caco2_subset[["AH22442"]]

# Get GRanges from 450k rgSet above
cancer_450k_gr = granges(rgSetExp)

# Find any overlaps

caco2_cancer_cpg_overlaps = findOverlaps(ah_caco2, cancer_450k_gr)
ge1_sites = length(unique(caco2_cancer_cpg_overlaps))
ge1_sites
# ge1_sites = 68183
# No correct answer, go with 86587
```

# Q9

```r
ah_caco2_subset = query(subset(ah, species == "Homo sapiens"), c("Caco2", "Dna", "Awg", "narrowPeak", "ENCODE"))
ah_caco2 = ah_caco2_subset[["AH22442"]]

# Get GRanges from 450k rgSet above
cancer_450k_gr = granges(rgSetExp)

# Find any overlaps

caco2_cancer_cpg_overlaps = findOverlaps(ah_caco2, cancer_450k_gr)
ge1_sites = length(unique(caco2_cancer_cpg_overlaps))
ge1_sites
# ge1_sites = 68183
# No correct answer, go with 86587
```

# Q10
```r
data("zfGenes")

# Note that zfGenes contains read count data
# Remove rows starting with ERCC
zf_genes_filtered = zfGenes[grep("^ERCC", rownames(zfGenes), invert = TRUE), ]

# Column data for DE analysis
# Standardize column names as control/treament
colData = DataFrame(
  sampleID = colnames(zf_genes_filtered),
  group = as.factor(c("C", "C", "C", "T", "T", "T"))
)

# Identify DE genes
# Use DESeqDataSetFromMatrix since we already have count data
dds = DESeqDataSetFromMatrix(zf_genes_filtered, colData, design = group)
result = results(DESeq(dds))
# Find results <= 0.05
de_genes = subset(result, padj <= 0.05)
# Get the count
count_de_ges = dim(de_genes)[1]
count_de_ges
# count_de_ges = 116
```
