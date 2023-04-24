using Plots
using CSV
using DataFrames
using LinearAlgebra
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

for file in glob("two_body/*.dat")
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


plot(x1[:, 1], x2[:, 1], label="m1")
plot!(x1[:, 2], x2[:, 2], label="m2")
xlabel!("x1 (pc)")
ylabel!("y1 (pc)")

savefig("two_body_xy.pdf")
