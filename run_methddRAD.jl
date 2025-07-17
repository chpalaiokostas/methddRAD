using Pkg
Pkg.activate(".") 

using methddRAD

merged_bam =ARGS[1] # merged bam file

bam_dir = ARGS[2] # directory of bam files
if isdir(bam_dir)
    bams = readdir(bam_dir)
else
    println("Need the directory with the bam files")
    exit()
end

samples = readlines(ARGS[3]) # file with sample names

merged_catalog(merged_bam)
catalog = "catalog_genomic_locations.bed" # expected name from merged_catalog

try 
    isfile(catalog)
    println("✅ Successfully created the catalog with the genomic locations")
    catalog_locations = DataFrame(CSV.File(catalog,header=[:Chrom,:Start,:End,:Strand,:Range]))
catch e
    println("❌ Error during catalog creation", e)
    exit(1)
end

# create pseudo bed files for each sample
open(bams, "r") do io
    bams = readlines(io)
    for bam in bams
        bam_to_sorted_bed(bam)
    end
end

matrix = zeros(Int,nrow(catalog_locations),length(samples)) 
df_counts = hcat(catalog_locations,DataFrame(matrix,samples))

#features = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
#                        for record in BED.Reader(open(catalog)))

features = Dict(key => 0 for key in df_counts[!,:Range])

for sample in samples
    feature_counts!(sample,catalog)
end

CSV.write("feature_counts_raw.txt",df_counts,delim="\t")

normalize_plus_bismark(df_counts)