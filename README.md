# Genomic Data Science Specialization

This repository contains the coursework for the Genomics data science course offered by John Hopkins University (through coursera).

Please see https://www.coursera.org/specializations/genomic-data-science for course details.

Certificate: https://www.coursera.org/account/accomplishments/specialization/KYQV4JH86EA8


## Setup Python Environment

### Create the conda environment
```shell
conda create -n genomics_data_science python=3.5
conda activate genomics_data_science
```

### Install the package
```shell
pip install -e .
```

### Additional dependencies


```shell
conda install -c bioconda samtools=1.2
conda install -c bioconda bedtools=2.25
conda install -c bioconda bowtie2=2.2.5
conda install -c bioconda bcftools=1.2
```

## Setup R Environment
All experiments were done using R studio using R version 4.1.3.
