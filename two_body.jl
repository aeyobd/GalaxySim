using Plots
using CSV
using DataFrames
using ProgressMeter
# using Printf
using LinearAlgebra
using StatsPlots

using GalaxySim

fpath = "../two_body/"
skip = 10

function get_col(col)
    file = fpath * col * ".dat"
    df = Array(CSV.read(file, DataFrame, header=false))[1:skip:end, 1:end-1]
    return df
end

x1s = get_col("x1")
x2s = get_col("x2")
x3s = get_col("x3")
v1s = get_col("v1")
v2s = get_col("v2")
v3s = get_col("v3")


ρs = get_col("rho")
Ts = log10.(get_col("T"))
ts = get_col("t")
du_C = get_col("du_C")
du_P = get_col("du_P")

dv_P = get_col("dv_P")
dv_G = get_col("dv_G")

m_star = get_col("mstar")

hs = get_col("h")
dts = get_col("dt")
Nn = get_col("N_neighbors")


Rs = @. sqrt(x1s^2 + x2s^2 + x3s^2)
Vs = @. sqrt(v1s^2 + v2s^2 + v3s^2);
du = du_C .+ du_P
dv = dv_P .+ dv_G;
energy = CSV.read("../two_body/energy.dat", DataFrame)[1:skip:end, :];
energy[!, "t"] = ts[:, 1];

plot(x1s[:, 1], x2s[:, 1])
plot!(x1s[:, 2], x2s[:, 2])


11.429

maximum(x2s[:, 1])*2

plot(v1s[:, 1], v2s[:, 1])
plot!(v1s[:, 2], v2s[:, 2])


length(ts)

energy

@df energy plot(:t, :tot, label="Total")
@df energy plot!(:t, :grav, label="gravitational")
@df energy plot!(:t, :kinetic, label="kinetic")
@df energy plot!(:t, :thermal, label="thermal")

histogram(log10.(dts[end, :]), range=(1, 10), bins=5)


