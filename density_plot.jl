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
    scatter!(df1.h, df1.f)
    hline!([0])
    xlims!(0, 40)
    ylims!(-0.01, 0.01)

    savefig("density.png")
end


make_plot()
