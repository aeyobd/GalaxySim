using Plots
using CSV
using DataFrames
using ProgressMeter
using Printf
using LinearAlgebra
using Logging
using LaTeXStrings


fpath = "./two_body_p/"
skip = 1

function get_col(col)
    file = fpath * col * ".dat"
    df = Array(CSV.read(file, DataFrame, header=false)[1:skip:end, 1:end-1])
    return df
end

x1s = get_col("x1")
v1s = get_col("x1")
ρs = get_col("rho")
Ts = (get_col("T"))
P = get_col("P")
ts = get_col("t")

t = (ts[:, 1] .+ ts[:, 2]) ./ 2 ./ 1e6

y = x1s
plot(t, y[:, 1], label="m1")
plot!(t, y[:, 2], label="m2")
xlabel!("time (Myr)")
ylabel!("position (pc)")
xlims!(0, 300)
ylims!(-50, 50)
savefig("two_body_p_pos.pdf")


y = ρs
plot(t, y[:, 1], label="m1")
plot!(t, y[:, 2], label="m2")
xlabel!("time (Myr)")
ylabel!(L"density (cm$^{-3}$)")
xlims!(0, 300)
savefig("two_body_p_rho.pdf")

y = Ts
plot(t, y[:, 1], label="m1")
plot!(t, y[:, 2], label="m2")
xlabel!("time (Myr)")
ylabel!("temperature (K)")
xlims!(0, 300)
savefig("two_body_p_T.pdf")
