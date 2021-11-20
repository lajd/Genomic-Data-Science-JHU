## Setup

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
