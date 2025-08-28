export count_meth_sites!

"""
    filter_and_count_records(bed_file::String, cut_sites::Dict{String, Set{Int}}, output_dir::String, counts::Dict)

Filters a pseudo BED file and simultaneously counts records associated with each cut site.
A record is counted for a cut site if its start or end position matches it.
If a record matches two sites (start and end), it is counted for both.
"""
function count_meth_sites!(sample::AbstractString, cut_sites::Set, df::DataFrame)
    sample_name = string(sample,".bed")
    counts = Dict()
    temp_dict = Dict(record => 0 for record in (cut_sites))
    reader = open(BED.Reader, sample_name)
    println("Processing $sample_name")
    for record in reader 
        is_match_found = false
        #start_site = string(BED.chrom(record),":",)
        chr = BED.chrom(record)
        start_site = BED.chromstart(record) + 1
        chr_start = string(chr,":",start_site)
        if chr_start in cut_sites
            temp_dict[chr_start] += 1
        end

        end_site = BED.chromend(record) - 1
        chr_end = string(chr,":",end_site)
        if chr_end in cut_sites
                temp_dict[chr_end] += 1
        end
    end
    close(reader)
    temp_df = DataFrame(Chrom=[split(key,":")[1] for key in keys(temp_dict)],
                        Pos=[split(key,":")[2] for key in keys(temp_dict)], 
                        Counts=[temp_dict[key] for key in keys(temp_dict)])
    sort!(temp_df,[:Chrom,:Pos])
    rename!(temp_df, "Counts" => sample)
    df[!,sample] = temp_df[!,sample]
end

