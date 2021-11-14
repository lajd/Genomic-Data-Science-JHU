## Setup

### Create the conda environment
```shell
conda create -n genomics_data_science python=3.9
conda activate genomics_data_science
```

### Install the package
```shell
pip install -e .
```

### Additional dependencies


```shell
conda install -c bioconda samtools
conda install -c bioconda bedtools
```