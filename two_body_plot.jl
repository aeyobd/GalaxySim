using Plots
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots


fpath = "two_body/"
skip = 10

function get_col(col)
    file = fpath * col * ".dat"
    df = Array(CSV.read(file, DataFrame, header=false, skipto=2))[1:skip:end, 1:end-1]
    return df
end

x1s = get_col("x1")
x2s = get_col("x2")
v1s = get_col("v1")
v2s = get_col("v2")

ts = get_col("t")

plot(x1s[:, 1], x2s[:, 1], label="m1")
plot!(x1s[:, 2], x2s[:, 2], label="m2")
xlabel!("x1 (pc)")
ylabel!("y1 (pc)")

savefig("two_body_xy.pdf")
