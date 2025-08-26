module methddRAD

using BED
using BioSequences
using CSV 
using DataFrames
using FASTX
using GenomicFeatures
using XAM

include("bam_to_bed.jl")
include("catalog_genomic_locations.jl")
include("feature_counts.jl")
include("ref_cut_sites.jl")
include("ref_count_meth_sites.jl")
include("normalize_plus_bismark.jl")
include("to_bismark.jl")


end 