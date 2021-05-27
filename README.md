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
sudo apt install -qq -y make cmake clang
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

In folder `Contrasts` run `snakemake`:
```
cd Contrasts
snakemake
```

## 3. Add features or debug in the python scripts
You made modifications to one of the python script, a notebook, this README.md, or you added new features.
You wish this work benefits to all (futur) users of this repository?
Please, feel free to open a [pull-request](https://github.com/ThibaultLatrille/AdaptaPop/pulls)

## TO DO:
- ~~Snakefile for generating gene specific SFS and run DFEM~~
- ~~Aggregate sites (SFS, fasta)~~
- ~~Histogram nearly-neutral (subsampling)~~
- ~~Histogram adaptive (boostrap)~~
- ~~Ontology test for genes~~
- ~~Control for omega in the nearly-neutral set (weighted sampling)~~
- ~~Detecting outliers based on CI~~*-
- ~~Regression on w_NA versus w_0 only for the nearly-neutral~~
- ~~Unfold SFS (needs polarization)~~
- ~~Get Fj-Fi for polymorphisms~~
- ~~Compress fasta results after DFEM runs~~
- Run on all samples
- Table with all statistics, each row is a population
- Nearly-neutral should filter adaptive sites in close proximity
- Run with w* (????)


From Martin’s talk in the morning we saw that positive linked selection influences piN/piS.
This means positive linked selection will influence the estimate of omega_NA in the MK setting.
Do you think that will bias your results?
Second/separate question, how do you think a fluctuation environment (fluctuating DFE) or recent change in population size will influence the congruence between the phylogenetic and the population genetic model?

How do you think the choice of the two different species is influencing the analysis ?
For example, for species within the same genus, or more distant.
How to choose that. And how does that affect?

Adaptation is actually computed using divergence in both approaches.
So it seems that what is congruent between phylogeny and polymorphism is w0 more than wA.
Am I right?

Assuming that everything hold and is robust, what would be you interpretation of the congruence?
Were you surprised by it?

About epistasis in the second part: could you select the genes detected as nearly neutral in the first part to get better estimates of Ne?
What you call epistatis is actually constraints on the evolution between sites, such as background selection would do, among other causes?

On the first part of your talk, could you use incongruence between omega_0 and omega_NA to learn about changes in Ne or changes in the DFE over time?

## Licence

The MIT License (MIT)

Copyright (c) 2019 Thibault Latrille

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


