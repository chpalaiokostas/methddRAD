export to_bismark

"""
    to_bismark(sample::AbstractString, df::DataFrame)
create bismark type of files. Expects a dataframe with normalized feature counts
"""
function to_bismark(sample::AbstractString, df::DataFrame)
    rel_meth = max_values .- df[!,sample]
    df = hcat(df[!,Cols(:Chr,:Start,:Strand)],rel_meth,df[!,sample],makeunique=true)
    df.pseudo_pattern_1 .= "CG"
    df.pseudo_pattern_2 .= "CG"
    CSV.write("bismark_files/bismark_format_$sample.txt", df; delim="\t", writeheader=false)
end