export count_meth_sites!

"""
    filter_and_count_records(bed_file::String, cut_sites::Dict{String, Set{Int}}, output_dir::String, counts::Dict)

Filters a pseudo BED file and simultaneously counts records associated with each cut site.
A record is counted for a cut site if its start or end position matches it.
If a record matches two sites (start and end), it is counted for both.
"""
function count_meth_sites!(sample::AbstractString, cut_sites::Dict{AbstractString, Vector{Int}}, df::DataFrame)
    sample_name = string(sample,".bed")
    counts = Dict()
    temp_dict = Dict(string(key,":",value) => 0 for (key,value) in cut_sites)
    #counts[sample_name] = Dict{Tuple{String, Int}, Int}()

    #output_filename = joinpath(output_dir, "filtered_" * sample_name)

    # BED.Writer ensures proper formatting of the output file
    #out_stream = open(BED.Writer, output_filename, "w")
    
    # BED.Reader provides an efficient, iterable stream of records
    open(BED.Reader, sample_name) do sample_features
        is_match_found = false
        chr = BED.chrom(sample_features)

        if haskey(cut_sites, chr)
            chr_sites = cut_sites[chr]

            start_site = BED.chromstart(sample_features) + 1
            if start_site in chr_sites
                site_key = string(chr,":", start_site)
                temp_dict[site_key] += 1
            end

            end_site = BED.chromend(record) - 1
            if end_site in chr_sites
                site_key = (chr, end_site)
                temp_dict[site_key] += 1
            end
        end
    end
    temp_df = DataFrame(Chrom=[split(key,":")[1] for key in keys(temp_dict)],
                        Pos=[split(key,":")[2] for key in keys(temp_dict)], 
                        Counts=[temp_dict[key] for key in keys(temp_dict)])
    sort!(temp_df,[:Chrom,:Pos])
    rename!(temp_df, "Counts" => sample)
    df[!,sample] = temp_df[!,sample]
end

