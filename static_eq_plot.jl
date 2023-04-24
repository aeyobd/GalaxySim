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
# Additionally, energy is mostly conserved.


using Plots
using CSV
using DataFrames
using LaTeXStrings
using StatsPlots
using Glob

skip = 10
function get_col(file)
    df = Array(CSV.read(file, DataFrame, 
            header=false, 
            comment="#", 
            delim=' ', 
            ignorerepeated=true)
        )[1:skip:end, :]
    return df
end

# a fancy loop to read in all the files
for file in glob("./static_eq/*.dat")
    fname, _ = splitext(basename(file))
    if fname == "energy"
        global energy
        energy = CSV.read(file, DataFrame, 
            header=["t", "thermal", "kinetic", "grav", "total"], 
            comment="#", 
            delim=' ', 
            ignorerepeated=true
            )[1:skip:end, :]

    else
        var = Symbol(fname)
        @eval $var = get_col($file)
    end
end

R = @. sqrt(x1^2 + x2^2 + x3^2)
V = @. sqrt(v1^2 + v2^2 + v3^2);
du = du_C .+ du_P .+ du_visc
dv = dv_P .+ dv_G .+ dv_visc;
ts = energy.t/1e6;

Nt, N = size(x1)


i = 2
scatter(R[i, :], log10.(rho[i,:]), label="t=0")
ylims!(minimum(log10.(rho[i,:])), maximum(log10.(rho[i, :])))

i=Nt
scatter!(R[i, :], log10.(rho[i,:]), label="t=100Myr")


x_model = LinRange(0, 400., 400)
y_model = log10.(123 * x_model .^ -2)

plot!(x_model, y_model, label=L"$1/r^2$", order=2, lw=5, c="black")
ylabel!(L"$\log \rho / {\rm cm}^{-3}$")
xlabel!(L"$R$ (pc)")

savefig("static_eq_rho.pdf")



i = 2

scatter(R[i, :], V[i,:], label="t=0")
i=Nt

scatter!(R[i, :], V[i,:], label="t=100Myr")

ylabel!("v (km/s)")
xlabel!("R (pc)")

savefig("static_eq_v.pdf")



