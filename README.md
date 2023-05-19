**Genes and sites under adaptation at the phylogenetic scale also exhibit adaptation at the population-genetic scale**\
Thibault Latrille, Nicolas Rodrigue, Nicolas Lartillot,\
_Proceedings of the National Academy of Sciences_,
Volume 120, Issue 11, March 2023, Pages e2214977120,\
doi: [10.1073/pnas.2214977120](https://doi.org/10.1073/pnas.2214977120)

**Compiled binaries and instructions for BayesCode are available at [ThibaultLatrille/bayescode](https://github.com/ThibaultLatrille/bayescode)**

# AdaptaPop

This repository is meant to provide the necessary scripts and instructions to reproduce the figures shown in the manuscript.
The experiments can either run on a local computer or in a cluster configuration (slurm).
The experiments are meant to run on Linux/Unix/MacOS operating systems.

The pipeline consist of three main steps, each with its own folder and Snakemake file:
- [I. Polymorphism](https://github.com/ThibaultLatrille/AdaptaPop#i-polymorphism---download-and-filter-vcf-files)
- [II. Divergence](https://github.com/ThibaultLatrille/AdaptaPop#ii-divergence---run-bayescode-on-orthomam)
- [III. Contrast polymorphism and divergence](https://github.com/ThibaultLatrille/AdaptaPop#iii-run-global-analysis-contrasting-polymorphism-and-divergence)

Section I and II are independent of each other and can be run in parallel, while section III depends on the results of I and II.
Alternatively, the results of each section can be downloaded at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7543458.svg)](https://doi.org/10.5281/zenodo.7543458)

If problems and/or questions are encountered, feel free to [open issues](https://github.com/ThibaultLatrille/AdaptaPop/issues).

## 0. Local copy
Clone the repository and `cd` to the dir.
```
git clone https://github.com/ThibaultLatrille/AdaptaPop
cd AdaptaPop
```

## 1. Installation

### General dependencies

Install python3 (>=3.9) and python3 packages
```
sudo apt install -qq -y python3-dev python3-pip
pip3 install --user snakemake numpy==1.23 scipy matplotlib pandas ete3 bio statsmodels seaborn rpy2 
```

### I. Polymorphism - Download and filter .vcf files 

Install bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
```
sudo apt install bedtools
```

In folder `Polymorphism` run `snakemake` in each subfolder:
```
cd Polymorphism
for FOLDER in ./*/ do 
    cd $FOLDER
    snakemake -j 8
    cd ..
done
```

Alternatively, you can download the filtered .vcf files at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7543458.svg)](https://doi.org/10.5281/zenodo.7543458) (`Polymorphism.zip`) and unzip them in the folder `Polymorphism`.

### II. Divergence - Run BayesCode on OrthoMam
This section is independent of the previous one ([I. Polymorphism](https://github.com/ThibaultLatrille/AdaptaPop#i-polymorphism---download-and-filter-vcf-files)).

To run on OrthoMam (mammalian orthologs), you must download the alignments and trees at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7543458.svg)](https://doi.org/10.5281/zenodo.7543458) (`OrthoMam.zip`) and unzip them in the folder `OrthoMam`. The input data will be already filtered and ready to be used in the folder `OrthoMam/Datasets`.

Install *BayesCode* from https://github.com/ThibaultLatrille/bayescode (see instructions there).

```
cd Orthomam
```

Run using snakemake, this requires access to large computation facilities:
```
snakemake -j 128
```
Or use the script `snakeslurm.sh` if run on a cluster (slurm) to submit the jobs.

Alternatively, you can download the results of BayesCode on OrthoMam at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7543458.svg)](https://doi.org/10.5281/zenodo.7543458) (`OrthoMam.zip`) and unzip them in the folder `OrthoMam`. The results will be located in the folder `OrthoMam/Experiments`.
Moreover, you can use this pipeline as a stand-alone to run `BayesCode` on your own set of alignments and trees by changing the `Snakefile` and/or the `config.yaml`.

### III. Run global analysis contrasting polymorphism and divergence

This section depends on the results of [I. Polymorphism](https://github.com/ThibaultLatrille/AdaptaPop#i-polymorphism---download-and-filter-vcf-files) and [II. Divergence](https://github.com/ThibaultLatrille/AdaptaPop#ii-divergence---run-bayescode-on-orthomam).

Install PAML
```
sudo apt install paml
```

In folder `Contrasts` run `snakemake`, this requires access to large computation facilities:
```
cd Contrasts
snakemake -j 128
```
Or use the script `snakeslurm.sh` if run on a cluster (slurm) to submit the jobs.

Alternatively, at the gene level, the rate of adaptation at the phylogenetic scale and at the population scale (McDonald & Kreitman) can be downloaded at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7543458.svg)](https://doi.org/10.5281/zenodo.7543458) (`GeneTable.tsv` and `MK_statistics.gz`).

## 3. Add features or debug in the python scripts
You made modifications to one of the python script, a notebook, this README.md, or you added new features.
You wish this work benefits to all (futur) users of this repository?
Please, feel free to open a [pull-request](https://github.com/ThibaultLatrille/AdaptaPop/pulls)

## Licence

The MIT License (MIT)

Copyright (c) 2019 Thibault Latrille

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


