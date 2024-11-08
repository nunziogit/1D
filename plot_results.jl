# plot_results.jl
function plot_results(xc, h_i, h, u_i, u, lx, timeout)
    p = plot(layout = (2, 1), size = (1200, 800))
    plot!(p[1], xc, h_i; linecolor = :black, linewidth = 2, label = "h_i",
        xlims = (0, lx), ylims = (minimum(h) - 0.1, maximum(h) + 0.1), xlabel = "x", ylabel = "h")
    plot!(p[1], xc, h; seriestype = :scatter, markercolor = :red, markersize = 5, label = "h",
        title = "time = $(round(timeout, digits = 1))")
    plot!(p[2], xc, u_i; linecolor = :black, linewidth = 2, label = "u_i",
        xlims = (0, lx), ylims = (minimum(u) - 0.1, maximum(u) + 0.1),
        xlabel = "x", ylabel = "u")
    plot!(p[2], xc, u; seriestype = :scatter, markercolor = :red, markersize = 5, label = "u")
    display(p)
end
