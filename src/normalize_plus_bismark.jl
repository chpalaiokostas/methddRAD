using CSV 
using DataFrames

include("to_bismark.jl")

feature_counts = CSV.read("feature_counts.txt", DataFrame)
sample_names = names(feature_counts, Not(r"Range|Chrom|Start|End"))

meth_counts_norm = mapcols(col -> col .* 1_000_000/sum(col), feature_counts, cols=sample_names)
meth_counts_norm_filt = meth_counts_norm[sum.(eachrow(meth_counts_norm[!, sample_names])) .> 0, :]
mapcols!(col -> Int.(round.(col;digits=0)), meth_counts_norm_filt, cols=sample_names)
max_values = maximum.(eachrow(meth_counts_norm_filt[!, sample_names]))

for sample in sample_names
    to_bismark(sample,meth_counts_norm_filt)
end