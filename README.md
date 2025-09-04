## methddRAD

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### Description

methddRAD has been designed for the analysis of reduced representation sequencing data that were created using methylation sensitive restriction enzymes. Those type of data offer a cost-effective approach for getting methylation related information in addition to the typical applications of reduced representation sequencing data. For the latter applications the interested user is advised to use a software like [Stacks](https://catchenlab.life.illinois.edu/stacks/)

The current version of methddRAD has the capacity of only handling data from reduced representation sequencing libraries that were created with the AciI restriction enzyme (recognition pattern: C|CGC).

Additional information about the library construction for the creation of data suitable to be analyzed with methddRAD can be found at: DOI: dx.doi.org/10.17504/protocols.io.kxygx43bol8j/v1 (https://www.protocols.io/view/a-cost-effective-protocol-for-methylome-and-genome-g73ubzqnx.html)


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
```
