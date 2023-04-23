# static_eq_plot
#
# A plot of the static equilibrium data
#
# Author Daniel Boyea (boyea.2@osu.edu)
# Created 21-04-2023
#
# As expected, the velocity is approximantly constant and 
# density is maintained as a 1/r^2 law.
# This means the density esimation is not perfect but qualitatively
# (<10%) accurate.
#


using Plots
using CSV
using DataFrames
using LaTeXStrings
using StatsPlots

fpath = "static_eq/"
skip = 10

function get_col(col)
    file = fpath * col * ".dat"
    df = Array(CSV.read(file, DataFrame, header=false, skipto=2))[1:skip:end, 1:end-1]
    return df
end


x1s = get_col("x1")
x2s = get_col("x2")
x3s = get_col("x3")
v1s = get_col("v1")
v2s = get_col("v2")
v3s = get_col("v3")



ρs = log10.(get_col("rho"))
Ts = log10.(get_col("T"))
ts = get_col("t")

hs = get_col("h")
dts = get_col("dt")
Nn = get_col("N_neighbors")

du_P = get_col("du_P")
du_C = get_col("du_C")
du_visc = get_col("du_visc")
du = @. du_P + du_C + du_visc


dv_P = get_col("dv_P")
dv_G = get_col("dv_G")
dv_visc = get_col("dv_visc")
dv = @. dv_P + dv_G + dv_visc


Rs = @. sqrt(x1s^2 + x2s^2 + x3s^2)
Vs = @. sqrt(v1s^2 + v2s^2 + v3s^2);
energy = CSV.read("static_eq/energy.dat", DataFrame)[1:skip:end, :];
energy[!, "t"] = ts[:, 1];

Nt, N = size(x1s)


i = 2
scatter(vec(Rs[i, :]), vec(ρs[i,:]), label="t=0")
ylims!(minimum(ρs[i,:]), maximum(ρs[i, :]))

i=Nt

scatter!(vec(Rs[i, :]), vec(ρs[i,:]), label="t=100Myr")


x_model = LinRange(0, 200., 100)
y_model = log10.(400 * x_model .^ -2)

plot!(x_model, y_model, label=L"$1/r^2$", order=2, lw=5, c="black")
ylabel!(L"$\log \rho / {\rm cm}^{-3}$")
xlabel!(L"$R$ (pc)")

savefig("static_eq_rho.pdf")



i = 2

scatter(vec(Rs[i, :]), vec(Vs[i,:]), label="t=0")
i=Nt

scatter!(vec(Rs[i, :]), vec(Vs[i,:]), label="t=100Myr")

ylabel!("v (km/s)")
xlabel!("R (pc)")

savefig("static_eq_v.pdf")



