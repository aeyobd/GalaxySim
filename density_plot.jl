# density_plot.jl
#
# Creates the plot of the density estimation
#
# Created 21-04-2023
# Author Daniel Boyea (boyea.2@osu.edu)
#
# As expected, the step size decreases as
# we approach the zero on f
#
# It is also promising that the derivative df
# is always negative and the function f
# is monotonically decreasing in this case.
# Essentially, the method of solving for density
# is okay.



using Plots
using LaTeXStrings
using CSV


function make_plot()
    df = CSV.File("density_f.dat")
    plot(df.h, df.f, label="f")
    plot!(df.h, df.df, label="df")
    xlabel!(L"$h$ (pc)")
    ylabel!(L"$f(h) = \rho - \rho_{\rm new}$")
    ylims!(-1, 1)
    df1 = CSV.File("density.dat")
    scatter!(df1.h, df1.f, label="calculated root")
    hline!([0], label=L"$f(h)=0$")
    xlims!(20, 60)
    ylims!(-0.0005, 0.0005)

    savefig("density.pdf")
end


make_plot()
