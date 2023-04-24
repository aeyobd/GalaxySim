# two_body_p_plot
# 
# creates the plots from the two_body_p.jl simulation
#
# Created 21-04-2023 
# Author Daniel Boyea (boyea.2@osu.edu)
#
# Writes plots to two_body_p_pos.pdf
# and two_body_p_rho.pdf
# and two_body_p_T.pdf
#
# The idea is that as the two particles approach eachother,
# density and temperature increase as they slow down due to 
# pressure, then the force of pressure causes the particles to 
# stop than accelerate away from eachother.
#
# So, qualitatively, this is working. (Unfortunantly, I haven't done 
# the analytic calculations.)


using Plots
using CSV
using DataFrames
using LinearAlgebra
using Logging
using LaTeXStrings
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
for file in glob("two_body_p/*.dat")
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


ts = energy.t/1e6

y = x1
plot(ts, y[:, 1], label="m1")
plot!(ts, y[:, 2], label="m2")
xlabel!("time (Myr)")
ylabel!("position (pc)")
xlims!(0, 300)
ylims!(-50, 50)
savefig("two_body_p_pos.pdf")


y = rho
plot(ts, y[:, 1], label="m1")
plot!(ts, y[:, 2], label="m2")
xlabel!("time (Myr)")
ylabel!(L"density (cm$^{-3}$)")
xlims!(0, 300)
savefig("two_body_p_rho.pdf")

y = T
plot(ts, y[:, 1], label="m1")
plot!(ts, y[:, 2], label="m2")
xlabel!("time (Myr)")
ylabel!("temperature (K)")
xlims!(0, 300)
savefig("two_body_p_T.pdf")
