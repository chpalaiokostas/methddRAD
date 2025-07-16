module methddRAD

using BED
using CSV 
using DataFrames
using GenomicFeatures
using XAM

include("catalog_genomic_locations.jl")
include("bam_to_bed.jl")
include("feature_counts.jl")
include("normalize_plus_bismark.jl")
include("to_bismark.jl")

merged_bam =ARGS[1]
samples = readlines(ARGS[2])
#bam_list = ARGS[1]

raw_catalog_locations(merged_bam)
catalog = "catalog_features.bed"
catalog_locations = DataFrame(CSV.File(catalog,header=[:Chrom,:Start,:End]))
catalog_locations.Range = string.(catalog_locations.Chrom,":",catalog_locations.Start,"-",catalog_locations.End)
select!(catalog_locations,:Range,Not(:Range))
matrix = zeros(Int,nrow(catalog_locations),length(samples)) 
df_counts = hcat(catalog_locations,DataFrame(matrix,samples))

feature_counts = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
                        for record in BED.Reader(open("catalog_features.bed")))

# bam_to_bed.jl
open(bam_list, "r") do io
    bam_files = readlines(io)
    for bam in bam_files
        bam_to_sorted_bed(bam)
    end
end

for sample in samples
    feature_counts!(sample)
end

CSV.write("feature_counts.txt",df_counts,delim="\t")

feature_counts = CSV.read("feature_counts.txt", DataFrame)
sample_names = names(feature_counts, Not(r"Range|Chrom|Start|End"))


end # module methddRAD
