using Plots
using CSV


df = CSV.File("density_f.dat")

plot(df.h, df.f, label="f", xaxis=:log)
plot!(df.h, df.df, label="df", xaxis=:log)

