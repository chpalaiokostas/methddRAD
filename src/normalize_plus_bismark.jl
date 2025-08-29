export normalize_plus_bismark
export normalize_cut_sites_bismark

"""
    normalize_plus_bismark(df::DataFrame)
Normalize read counts based on number of reads each sample has.
The normalized value is multiplied by 1_000_000 and rounded to nearest integer.
"""
function normalize_plus_bismark(df::DataFrame)
    sample_names = names(df, Not(r"Chrom|End|Range|Start|Empty|Strand"))
    meth_counts_norm = mapcols(col -> col .* 1_000_000/sum(col), df, cols=sample_names)
    meth_counts_norm_filt = meth_counts_norm[sum.(eachrow(meth_counts_norm[!, sample_names])) .>= 10, :]
    mapcols!(col -> Int.(round.(col;digits=0)), meth_counts_norm_filt, cols=sample_names)
    max_values = maximum.(eachrow(meth_counts_norm_filt[!, sample_names]))
    for sample in sample_names
        to_bismark(sample,max_values,meth_counts_norm_filt)
    end
end

function normalize_cut_sites_bismark(df::DataFrame)
    sample_names = names(df, Not(r"Chrom|Pos"))
    df_filtered = df[sum.(eachrow(df[!, sample_names])) .>= 10, :]
    #meth_counts_norm_filt = meth_counts_norm[sum.(eachrow(meth_counts_norm[!, sample_names])) .>= 10, :]
    meth_counts_norm = mapcols(col -> col .* 1_000_000/sum(col), df_filtered, cols=sample_names)
    #meth_counts_norm_filt = meth_counts_norm[sum.(eachrow(meth_counts_norm[!, sample_names])) .>= 10, :]
    mapcols!(col -> Int.(round.(col;digits=0)), meth_counts_norm, cols=sample_names)
    max_values = maximum.(eachrow(meth_counts_norm[!, sample_names]))
    for sample in sample_names
        ref_cut_sites_to_bismark(sample,max_values,meth_counts_norm)
    end
end