using ArgParse
using Pkg
Pkg.activate(".") 

using methddRAD

function parse_commandline()
    s = ArgParseSettings(description="Get methylation estimates from methddRAD type of data")

    @add_arg_table! s begin
        "--merged_bam","-m"
            help = "Provide the path to the merged bam file."
            required = true
            arg_type = String
        "--bam_files","-b"
            help = "Provide the directory where the bam files are located."
            required = true
            arg_type = String
        "--samples","-s"
            help = "List of sample names"
            required = true
            arg_type = String
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    merged_bam = parsed_args["merged_bam"] # merged bam file

    bam_dir = parsed_args["bam_files"] # directory of bam files
    if isdir(bam_dir)
        bams = readdir(bam_dir)
        filter!(x -> endswith(x,"bam"),bams)
    else
        println("Need the directory with the bam files")
        exit()
    end

    samples = readlines(parsed_args["samples"]) # file with sample names

    merged_catalog(merged_bam)
    catalog = "catalog_genomic_locations.bed" # expected name from merged_catalog

    try 
        isfile(catalog)
        println("✅ Successfully created the catalog with the genomic locations")
        catalog_locations = DataFrame(CSV.File(catalog,header=[:Chrom,:Start,:End,:Strand,:Range]))
    catch e
        println("❌ Error during catalog creation", e)
        exit(1)
    end

    # create pseudo bed files for each sample
    open(bams, "r") do io
        bams = readlines(io)
        for bam in bams
            bam_to_sorted_bed(bam)
        end
    end

    matrix = zeros(Int,nrow(catalog_locations),length(samples)) 
    df_counts = hcat(catalog_locations,DataFrame(matrix,samples))

    #features = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
    #                        for record in BED.Reader(open(catalog)))

    for sample in samples
        feature_counts!(sample,catalog,df_counts)
    end

    CSV.write("feature_counts_raw.txt",df_counts,delim="\t")

    normalize_plus_bismark(df_counts)
end

main()