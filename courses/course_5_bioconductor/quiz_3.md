# Quiz 3

# Setup

```r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biocLite")
BiocManager::install("biomaRt")

# ALL experimental dataset
install("ALL")
library("ALL")
data(ALL)

install("biomaRt")
library("biomaRt")
# Load Microarray for ALL
install("hgu95av2.db")
library("hgu95av2.db")

install("airway")
library(airway)
data(airway)

library(AnnotationHub)
ah = AnnotationHub()


install("minfiData")
install("minfi")
library(minfiData)
library(minfi)

install("GEOquery")
library(GEOquery)

library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
```

# Q1

```r

mean_expr_col_5 = mean(exprs(ALL)[, 5])
mean_expr_col_5
# mean_expr_col_5 = 5.629627
```

# Q2
```r
mart = useMart(host="feb2014.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("hsapiens_gene_ensembl", mart)

# Get the ALL annotations
annotations = getBM(
     attributes=c("affy_hg_u95av2", "ensembl_gene_id"),
     filters="affy_hg_u95av2", values=featureNames(ALL),
     mart=ensembl
)
probsets_with_ge2_gene_ids = sum(table(annotations[, 1]) > 1)
# probsets_with_ge2_gene_ids = 1045

```

# Q3

```r

annotations_autosomes <- getBM(
    attributes=c("affy_hg_u95av2", "ensembl_gene_id", "chromosome_name"),
    filters=c("affy_hg_u95av2", "chromosome_name"), values=list(featureNames(ALL), c(1:22)),
    mart=ensembl
)

probsets_autosomes_with_ge2_gene_ids = sum(table(annotations_autosomes[,1])>1)
probsets_autosomes_with_ge2_gene_ids
# probsets_with_autosomes_ge2_gene_ids = 2234


probsets_autosomes_with_ge1_genes = length(unique(annotations_autosomes[, 2]))
probsets_autosomes_with_ge1_genes
# probsets_autosomes_with_ge1_genes = 11016

```

# Q4

```r
colnames(MsetEx)
# [1] "5723646052_R02C02" "5723646052_R04C01" "5723646052_R05C02" "5723646053_R04C02" "5723646053_R05C02"
# [6] "5723646053_R06C02"
# Desired column 5723646052_R04C01 is column # 1 (or 2 if including the index)
MsetEx["5723646052_R04C01]

# Get the data
methData = getMeth(MsetEx)
mean_meth_chl = mean(methData[, 2])
mean_meth_chl
# mean_meth_chl = 7228.277
```

# Q5

```r
gse = getGEO("GSE788")

gse_expr = exprs(gse[[1]])
colnames(gse_expr)
# [1] "GSM9011" "GSM9024" "GSM9025" "GSM9026" "GSM9027"
# [6] "GSM9028" "GSM9030" "GSM9031" "GSM9032" "GSM9033"
# GSM9024 is column 1 (col 2 including index)
mean_expr_gsm9024 = mean(gse_expr[,2])
mean_expr_gsm9024
# mean_expr_gsm9024 = 756.432
```

# Q6

```r
mean_avg_length = mean(airway$avgLength)
mean_avg_length
# mean_avg_length = 113.75
```

# Q7

```r
airway_assay = assay(airway)
# colnames(airway_assay)
# [1] "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513" "SRR1039516" "SRR1039517"
# [7] "SRR1039520" "SRR1039521"
# SRR1039512 is col 3
ensemble_gene_count_with_ge1_read = sum(airway_assay[, 3] > 0)
ensemble_gene_count_with_ge1_read
# ensemble_gene_count_with_ge1_read = 25699
```

# Q8

```r
hg19_exons = exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
standard_chromosome_exons = keepStandardChromosomes(hg19_exons, pruning.mode="coarse")
autosome_exons = dropSeqlevels(standard_chromosome_exons, c("chrX", "chrY", "chrM"), pruning.mode = "coarse")

# Need to map chromosome names from "chr<x>" -> "<x>", i.e. map to NCBI format
ncbi_chr_map = mapSeqlevels(seqlevels(autosome_exons), "NCBI")
autosome_exons_ncbi = renameSeqlevels(autosome_exons, ncbi_chr_map)

overlap_features_with_autosome_exons = subsetByOverlaps(airway, autosome_exons_ncbi)
dim(overlap_features_with_autosome_exons)[1]
# dim(overlap_features_with_autosome_exons)[1] = 26276
```

# Q9

```r

# colnames(airway)
# [1] "SRR1039508" "SRR1039509" "SRR1039512"
# [4] "SRR1039513" "SRR1039516" "SRR1039517"
# [7] "SRR1039520" "SRR1039521"

# SRR1039508 is col 1

airway_sample_srr1039508 = airway[, 1]
srr1039508_autosomes_ncbi_subset = subsetByOverlaps(airway_sample_srr1039508, autosome_exons_ncbi)

srr1039508_total_reads = sum(assay(airway_sample_srr1039508)[, 1])
srr1039508_overlap_reads = sum(assay(srr1039508_autosomes_ncbi_subset)[, 1])

frac_srr1039508_overlapping_autosome = srr1039508_overlap_reads/srr1039508_total_reads
frac_srr1039508_overlapping_autosome
# frac_srr1039508_overlapping_autosome = 0.9004193

```

# Q10

```r

ah_e096_record = query(ah, c("H3K4me3", "narrowPeak", "E096"))[["AH30596"]]
ah_e096_autosomes_ = keepStandardChromosomes(ah_e096_record, pruning.mode="coarse")
ah_e096_autosomes = dropSeqlevels(ah_e096_autosomes_, c("chrX", "chrY", "chrM"), pruning.mode = "coarse")
ah_e096_autosomes_ncbi = renameSeqlevels(ah_e096_autosomes, ncbi_chr_map)


ah_record_autosome = keepSeqlevels(ah_e096_record, autosome, pruning.mode = "coarse")
ah_record_ncbi = renameSeqlevels(ah_record_autosome, txdb_ncbi)


```

