using Pkg
Pkg.activate(".") 

using methddRAD

merged_bam =ARGS[1]
bams = readlines(ARGS[2])
samples = readlines(ARGS[3])

merged_catalog(merged_bam)
catalog = "catalog_genomic_locations.bed" # expected name from merged_catalog
catalog_locations = DataFrame(CSV.File(catalog,header=[:Chrom,:Start,:End,:Strand,:Range]))

# create pseudo bed files for each sample
open(bams, "r") do io
    bams = readlines(io)
    for bam in bams
        bam_to_sorted_bed(bam)
    end
end

matrix = zeros(Int,nrow(catalog_locations),length(samples)) 
df_counts = hcat(catalog_locations,DataFrame(matrix,samples))

features = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
                        for record in BED.Reader(open(catalog)))

for sample in samples
    feature_counts!(sample,catalog)
end

CSV.write("feature_counts_raw.txt",df_counts,delim="\t")

normalize_read_counts(df_counts)