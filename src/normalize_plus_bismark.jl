export normalize_plus_bismark

"""
    normalize_read_counts(df::DataFrame)
Normalize read counts based on number of reads each sample has.
The normalized value is multiplied by 1_000_000 and rounded to nearest integer.
"""
function normalize_plus_bismark(df::DataFrame)
    sample_names = names(feature_counts, Not(r"Chrom|End|Range|Start|Strand"))
    meth_counts_norm = mapcols(col -> col .* 1_000_000/sum(col), df, cols=sample_names)
    meth_counts_norm_filt = meth_counts_norm[sum.(eachrow(meth_counts_norm[!, sample_names])) .> 0, :]
    mapcols!(col -> Int.(round.(col;digits=0)), meth_counts_norm_filt, cols=sample_names)
    max_values = maximum.(eachrow(meth_counts_norm_filt[!, sample_names]))
    for sample in sample_names
        to_bismark(sample,meth_counts_norm_filt)
    end
end