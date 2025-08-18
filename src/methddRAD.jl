module methddRAD

#using BED
#using CSV 
#using DataFrames
#using GenomicFeatures
#using XAM

include("catalog_genomic_locations.jl")
include("bam_to_bed.jl")
include("feature_counts.jl")
include("normalize_plus_bismark.jl")
include("to_bismark.jl")

end # module methddRAD
