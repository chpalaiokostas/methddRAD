## methddRAD

### How to run the software

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
