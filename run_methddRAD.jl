using ArgParse

using Pkg
Pkg.activate(".")

using Base.Threads
using BED
using BioSequences
using CSV 
using DataFrames
using GenomicFeatures
using FASTX
using methddRAD
using XAM

"""
    parse_commandline()
Parses command-line arguments
"""
function parse_commandline()
    s = ArgParseSettings(description="Get methylation estimates from methddRAD type of data")

    @add_arg_table! s begin
        "--merged_bam","-m"
            help = "Provide the path to the merged bam file."
            required = false
            arg_type = String
        "--bam_files","-b"
            help = "Provide the directory where the bam files are located."
            required = true
            arg_type = String
        "--samples","-s"
            help = "List of sample names"
            required = true
            arg_type = String
        "--guided"
            help = "Prior information about the expected cutting sites in the genome"
            action = :store_true
        "--genome"
            help = "Provide the reference genome file"
            requred = false
            arg_type = String
    end
    return parse_args(s)
end

function main()
    println("✅ Arguments parsed. Activating environment and loading packages...")
    parsed_args = parse_commandline()

    if !parsed_args["guided"] && isnothing(parsed_args["merged_bam"])
        error("The --merged_bam is needed when not using --guided")
    end

    if parsed_args["guided"] && isnothing(parsed_args["genome"])
        error("The --genomed is needed when using --guided")
    end

    if !parsed_args["guided"]
        println("Will run a denovo approach to create a catalog with observed sites")
        merged_bam = parsed_args["merged_bam"] # merged bam file

        bam_dir = abspath(parsed_args["bam_files"]) # directory of bam files
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
        catalog_locations = nothing

        try 
            isfile(catalog)
            println("✅ Successfully created the catalog with the genomic locations")
            catalog_locations = DataFrame(CSV.File(catalog,header=[:Chrom,:Start,:End,:Range,:Empty,:Strand]))
        catch e
            println("❌ Error during catalog creation", e)
            exit(1)
        end

        # create pseudo bed files for each sample
        @threads for bam in string.(bam_dir,bams)
            bam_to_sorted_bed(bam)
        end

        matrix = zeros(Int,nrow(catalog_locations),length(samples)) 
        df_counts = hcat(catalog_locations,DataFrame(matrix,samples))
        features = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
                                for record in BED.Reader(open(catalog)))
        #features = Dict(key_name = string(seqname(record),":",leftposition(record),"-",rightposition(record)) => 0 
        #                    for record in BED.Reader(open("catalog_genomic_locations.bed")))
        for sample in samples
            feature_counts!(sample,catalog,df_counts,features)
        end

        CSV.write("feature_counts_raw.txt",df_counts,delim="\t")

        normalize_plus_bismark(df_counts)
    else
        fasta = parse_args["genome"]
        all_cut_sites(fasta)
        

    end
end

main()
