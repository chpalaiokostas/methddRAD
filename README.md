## methddRAD - A cost effective approach for deriving methylation related information

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### Description

methddRAD has been designed for the analysis of reduced representation sequencing data that were created using methylation sensitive restriction enzymes. Those type of data offer a cost-effective approach for getting methylation related information in addition to the typical applications of reduced representation sequencing data. For the latter applications the interested user is advised to use a software like [Stacks](https://catchenlab.life.illinois.edu/stacks/)

The current version of methddRAD has the capacity of only handling data from reduced representation sequencing libraries that were created with the AciI restriction enzyme (recognition pattern: C|CGC). Additional information about the library construction for the creation of data suitable to be analyzed with methddRAD can be found at: [DOI: dx.doi.org/10.17504/protocols.io.kxygx43bol8j/v1](https://www.protocols.io/view/a-cost-effective-protocol-for-methylome-and-genome-g73ubzqnx.html)

There are two approaches that can be used to analyzed the data:

**1st**: Here a similar approach to Marconi et al 2019 (https://www.nature.com/articles/s41598-019-51423-2) is followed based on creating a genome wide catalog of expected methylation sites where nearby sites are merged. Here the user need to create beforehand a merged bam file from all the available samples. Using samtools is probable the easiest way for this (https://www.htslib.org/doc/samtools-merge.html).

**2nd**: Here a reference genome will not to be supplied. methddRAD would scan the reference genome and identify the sites where the AciI enzyme will cut. Those sites will be used as the baseline for deriving methylation related information for each sample.

Uppon successful completion of a methddRAD run bismark type of files will be created for each sample containing information about the level of methylation of each site. Those files can be supplied to methylation dedicated software for downstream analysis.

### How to run the software

```
usage: julia ./run_methddRAD.jl [-m MERGED_BAM] -b BAM_FILES -s SAMPLES
                        [--guided] [--genome GENOME] [-h]

Get methylation estimates from methddRAD type of data

optional arguments:
  -m, --merged_bam MERGED_BAM
                        Provide the path to the merged bam file.
  -b, --bam_files BAM_FILES
                        Provide the directory where the bam files are
                        located.
  -s, --samples SAMPLES
                        List of sample names
  --guided              Prior information about the expected cutting
                        sites in the genome
  --genome GENOME       Provide the reference genome file
  -h, --help            show this help message and exit

To follow the 1st analysis approach approach use the flag *-m*. In this case the --guided and --genome flags are not to be used.

Likewise for the 2nd analysis approach one should not use the *-m* flag. Instead use the --guided and --genome ones.
```
