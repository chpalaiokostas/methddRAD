module methddRAD

include("catalog_genomic_locations.jl")
include("bam_to_bed.jl")
include("feature_counts.jl")
include("normalize_plus_bismark.jl")
include("to_bismark.jl")

merged_bam =ARGS[1]
samples = readlines(ARGS[2])
#bam_list = ARGS[1]

raw_catalog_locations(merged_bam)

# bam_to_bed.jl
open(bam_list, "r") do io
    bam_files = readlines(io)
    for bam in bam_files
        bam_to_sorted_bed(bam)
    end
end


end # module methddRAD
