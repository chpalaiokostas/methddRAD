export get_strand
export raw_catalog_locations
export merged_catalog

# didn't find any built-in function for finding the strand in XAM
"""
    get_strand(record::BAM.Record)
infers the strand of the sequenced record 
"""    
function get_strand(record::BAM.Record)
    flag = BAM.flags(record)
    if (flag & 0x10) != 0
        return "-" 
    else
        return "+" 
    end
end

"""
    raw_catalog_locations(reader::BAM.Reader)
Find all unique locations.
Needs a merged bam file from all available samples
"""
function raw_catalog_locations(reader::AbstractString)
    reader = open(BAM.Reader,reader)
    record = BAM.Record()
    annotation = DataFrame()
    # populate annotation with Chrom, Start, End, Strand
    while !eof(reader)
        empty!(record)
        read!(reader, record)
        if !BAM.ismapped(record)
            continue
        end
        if (BAM.refname(record) == BAM.nextrefname(record)) && (200 <= BAM.templength(record) <= 600)
            push!(annotation,(Chrom=BAM.refname(record), 
            Start=BAM.position(record) - 1,
            End=BAM.position(record) + BAM.templength(record) - 1,
            Strand=get_strand(record)))
        end
    end
    close(reader)
    unique!(annotation)
    CSV.write("genomic_locations_raw.csv", annotation)
end

"""
    merged_catalog(reader::BAM.Reader)
Creates pseudo bed file merging overlaps.
"""
function merged_catalog(reader::AbstractString)
    raw_catalog_locations(reader)
    annotation = CSV.read("genomic_locations_raw.csv", DataFrame)
    annotation_overlap = DataFrame()
    previous = nothing
    for current in eachrow(annotation)
        if !isnothing(previous) && current.Chrom == previous.Chrom && current.Start <= previous.End
            if (max(current.End, previous.End) - min(current.Start, previous.Start)) <= 800
                current.End = max(current.End, previous.End)
                current.Start = min(current.Start, previous.Start)
            else
                push!(annotation_overlap,previous)
            end    
        elseif !isnothing(previous) && current.Chrom == previous.Chrom && current.Start > previous.End
            push!(annotation_overlap,previous)
        end
        previous = current
    end
    push!(annotation_overlap,previous)
    annotation_overlap.Range = string.(annotation_overlap.Chrom,":",annotation_overlap.Start,
                                    "-",annotation_overlap.End)
    annotation_overlap.Empty .= 0
    select!(annotation_overlap,Not(:Strand),:Strand)
    sort!(annotation_overlap,[:Chrom, :Start])
    CSV.write("catalog_genomic_locations.bed", delim="\t", header=false,annotation_overlap)
end