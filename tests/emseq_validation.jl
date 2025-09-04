ENV["GKSwstype"] = "100"

using CSV, DataFrames
using Plots
using Statistics

emseq_samples = readlines("emseq_samples.txt")
rx = r"\w+-\d{4}-s-2020-(\w+?)_\w+"

df_summary = DataFrame()

for file in emseq_samples
    m = match(rx, file)
    meth_ddRAD_sample = string("bismark_format_em_seq_",m[1],".txt")
    meth_ddRAD = CSV.read(meth_ddRAD_sample,header=false, DataFrame)
    meth_ddRAD.Key = string.(meth_ddRAD.Column1,":",meth_ddRAD.Column2)
    emseq = CSV.read(file,header=false,DataFrame)
    emseq.Key = string.(emseq.Column1,":",emseq.Column2)
    to_validate = innerjoin(meth_ddRAD, emseq; on=[:Key],renamecols = "_left" => "_right")
    to_validate.emseq_reads = to_validate.Column5_right .+ to_validate.Column6_right
    to_validate = to_validate[to_validate.emseq_reads .>=10,:]
    to_validate.meth_ddRAD_perc = to_validate.Column4_left ./(to_validate.Column4_left .+ to_validate.Column5_left) .* 100
    to_validate.diff_emseq_methddRAD = abs.(to_validate.meth_ddRAD_perc .- to_validate.Column4_right)
    push!(df_summary,(Id=m[1],Meth_diff_median=median(to_validate.diff_emseq_methddRAD),N_sites=nrow(to_validate)))
    histogram(to_validate.diff_emseq_methddRAD; labels=false)
    savefig(string("meth_diff_",m[1],".pdf"))
end
    
CSV.write("Methylation_diff.txt",df_summary,delim="\t")