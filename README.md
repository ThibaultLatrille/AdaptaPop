# AdaptaPop

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.
The experiments can either run on a local computer or in a cluster configuration (slurm).

The experiments are meant to run on Linux/Unix/MacOS operating systems.

If problems and/or questions are encountered, feel free to [open issues](https://github.com/ThibaultLatrille/AdaptaPop/issues).

## 0. Local copy
Clone the repository and `cd` to the dir.
```
git clone https://github.com/ThibaultLatrille/AdaptaPop
cd AdaptaPop
```

## 1. Installation

### General dependencies

Install python3 packages
```
sudo apt install -qq -y python3-dev python3-pip
pip3 install snakemake scipy numpy matplotlib pandas ete3 bio statsmodels --user
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

### II. Divergence - Run BayesCode on OrthoMam

```
cd Orthomam
```

Install the compiling toolchains:
```
sudo apt install -qq -y make cmake clang openmpi-bin openmpi-common libopenmpi-dev
```
Clone and compile the C++ code for *BayesCode*
```
git clone https://github.com/ThibaultLatrille/bayescode && cd bayescode && git checkout dev && make release && cd ..
```
Run using snakemake, this requires access to large computation facilities:
```
snakemake -j 128
```
Or use the script `snakeslurm.sh` if run on a cluster (slurm) to submit the jobs.
### III. Run global analysis contrasting polymorphism and divergence

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


