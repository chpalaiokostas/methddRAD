using BED
using CSV
using DataFrames
using GenomicFeatures
using XAM

samples = readlines(ARGS[1])
catalog = "catalog_features.bed"
catalog_locations = DataFrame(CSV.File(catalog,header=[:Chrom,:Start,:End]))
catalog_locations.Range = string.(catalog_locations.Chrom,":",catalog_locations.Start,"-",catalog_locations.End)
select!(catalog_locations,:Range,Not(:Range))
matrix = zeros(Int,nrow(catalog_locations),length(samples)) 
df_counts = hcat(catalog_locations,DataFrame(matrix,samples))

feature_counts = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
                        for record in BED.Reader(open("catalog_features.bed")))

"""
    feature_counts!(sample::AbstractString)
find number of reads each sample has on the appropriate locations of the df_counts                     
"""
function feature_counts!(sample::AbstractString)
    temp_dict = copy(feature_counts)
    sample_name = string(sample,".bed")
    open(BED.Reader, sample_name) do features_x
        open(BED.Reader, "catalog_features.bed") do features_y
            for (x, y) in eachoverlap(features_x, features_y)
                key_name = string(seqname(y),":",leftposition(y),"-",rightposition(y))
                temp_dict[key_name] = get(temp_dict, key_name, 0) + 1    
            end
        end
    end
    temp_df = DataFrame(Range=[key for key in keys(temp_dict)],
                                Counts=[temp_dict[key] for key in keys(temp_dict)])
    transform!(temp_df, 
                :Range => ByRow(x -> split(x,":")[1]) => :Chrom,
                :Range => ByRow(x -> split(split(x,":")[2],"-")[1]) => :Start,
                :Range => ByRow(x -> split(split(x,":")[2], "-")[2]) => :End)
    temp_df.Start = parse.(Int, temp_df.Start)
    temp_df.End = parse.(Int, temp_df.End)
    sort!(temp_df,[:Chrom,:Start])
    rename!(temp_df, "Counts" => sample)
    df_counts[!,sample] = temp_df[!,sample]
end

for sample in samples
    feature_counts!(sample)
end

CSV.write("feature_counts.txt",df_counts,delim="\t")