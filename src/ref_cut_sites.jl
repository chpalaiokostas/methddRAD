export find_cut_sites
export all_cut_sites

"""
    find_cut_sites(sequence::LongDNASeq, recognition_site::ExactSearchQuery)

Finds all start positions of a recognition site within a larger DNA sequence.

# Arguments
- `sequence::LongDNASeq`: A sequence record from a FASTA file.
- `recognition_site::ExactSearchQuery`: The DNA sequence of the enzyme's recognition site.

# Returns
- `Vector{Int}`: A vector of 1-based start positions of all matches.
"""
function find_cut_sites(sequence::LongDNA{4}, recognition_site::ExactSearchQuery)
    positions = Int[]
    current_pos = 1
    while true
        match_range = findnext(recognition_site, sequence, current_pos)
        if match_range === nothing
            break
        else
            cut_site = first(match_range) + 1
            push!(positions, cut_site)
            current_pos = cut_site
        end
    end
    return positions
end

"""
    all_cut_sites(genome::AbstractString)

Finds all start positions of a recognition site within a genome.

# Arguments
- `genome::AbstractString`: The genome in the form of a FASTA file.
"""
function all_cut_sites(genome::AbstractString)
    enzyme_name = "AciI"
    recognition_site = dna"CCGC"
    println("Will check for cut sites in the reference genome....")
    println("Enzyme: $enzyme")
    println("Recognition Site: $recognition_site")
    println("Genome: $fasta\n")

    rev_comp_site = reverse_complement(recognition_site)
    is_palindromic = (recognition_site == rev_comp_site)

    if !is_palindromic
        println("Note: The site is not palindromic; searching both strands separately.")
    end

    all_sites_by_sequence = Dict{String, Vector{Int}}()
    total_sites_found = 0
    FASTAReader(open(genome)) do reader
        for record in reader
            genome_id = identifier(record)
            genome_seq = LongDNA{4}(sequence(record))

            # Find hits on the forward strand
            forward_hits = find_cut_sites(genome_seq, ExactSearchQuery(recognition_site))

            # Find hits on the reverse strand (unless palindromic)
            reverse_hits = Int[]
            if !is_palindromic
                reverse_hits = find_cut_sites(genome_seq, ExactSearchQuery(rev_comp_site))
            end

            # Combine, remove duplicates, and sort the positions for the current sequence
            combined_sites = sort(collect(Set([forward_hits..., reverse_hits...])))

            # If sites were found, store them in the dictionary
            if !isempty(combined_sites)
                all_sites_by_sequence[genome_id] = combined_sites
                total_sites_found += length(combined_sites)
            end
        end
    end
    if isempty(all_sites_by_sequence)
        println("\n❌ No cut sites found for $enzyme in any sequence.")
    else
        #num_sequences_with_hits = length(values(all_sites_by_sequence))
        println("\n✅ Analysis Complete. Found $total_sites_found total site(s).")
            
        cut_sites = DataFrame()
        for id in keys(all_sites_by_sequence)
            sites = all_sites_by_sequence[id]
            for site in sites
                push!(cut_sites, (Chrom=id, Pos=site))
            end            
        end
        sort!(cut_sites,[:Chrom,:Pos])
        CSV.write("cut_sites.csv", cut_sites)
    end
end