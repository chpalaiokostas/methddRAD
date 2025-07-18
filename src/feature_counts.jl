export feature_counts!

"""
    feature_counts!(sample::AbstractString)
find number of reads each sample has on the appropriate locations of the df_counts                     
"""
function feature_counts!(sample::AbstractString, catalog::AbstractString, df::DataFrame)
    #temp_dict = Dict(key => 0 for key in df[!,:Range])
    temp_dict = copy(features)
    sample_name = string(sample,".bed")
    open(BED.Reader, sample_name) do features_x
        open(BED.Reader, catalog) do features_y
            for (x, y) in eachoverlap(features_x, features_y)
                key_name = string(seqname(y),":",leftposition(y),"-",rightposition(y))
                temp_dict[key_name] += 1
                #temp_dict[key_name] = get(temp_dict, key_name, 0) + 1    
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
    df[!,sample] = temp_df[!,sample]
end

