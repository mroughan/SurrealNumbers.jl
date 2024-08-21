# do some plots from the output of random_surreals_sim.jl
#   this case does plots for the final distributions
#   in particular, this compares the Uniform and Binary splitting function
using DataFrames
using CSV
using Printf

using Random
using Distributions

using Plots
using PyPlot
pyplot()

input_dir = "Data/Compiled/"
plot_dir = "Plots/"

include("random_surreals.jl")
include("utilities.jl")
include("parameters.jl") 
include("plot_parameters.jl") 

α = 0.4; λ = 1.5
α = 0.8; λ = 3.5

g_0 = 1
n = 4000
Nseeds = 30
PopulationSize = n * Nseeds

distributions = "P-U"
base_file_name = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d", distributions, λ,  α , n, g_0, PopulationSize))

# plot final distributions
p1 = Dict{String,Any}() # scatter plot of the PMF
p2 = Dict{String,Any}() # log-y PMF
p3 = Dict{String,Any}() # log-log CCDF

x_thresh = 30

for k in ["Generation","Nodes","Edges","Log2_Denominator"]
    global input_file = "$(base_file_name)_$k.csv" 
    global DF = CSV.read(input_file, DataFrame)
    global x = DF[:, Symbol(k)]
    global y = DF[:, Symbol("f($k)")] ./ DF[:, Symbol("Population Size") ]

    if k=="Log2_Denominator"
        global label1 = "Denominator-exponent"
    else
        global label1 = "$k"
    end
    
    PyPlot.figure(k)
    if maximum(x) <= x_thresh
        global p1[k] = Plots.plot(x, y; seriestype=:scatter, 
                                    marker = marker[1], markercolor = mc[1], markerstrokecolor = mc[1],
                                    alpha=0.8,
                                    label="$label1 (Uniform)"
                                )
    else
        global p1[k] = Plots.plot(x, y; color = mc[1], label="$k (Uniform)")
    end
    Plots.plot!(p1[k]; xlabel="k", ylabel="PMF")
    if k=="Log2_Denominator"
        Plots.plot!(p1[k]; ylims=(0.0, 1.0) )
    end
    
    global distributions2 = "P-B"
    global base_file_name2 = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d", distributions2, λ,  α , n, g_0, PopulationSize))
    global input_file2 = "$(base_file_name2)_$k.csv" 
    global DF2 = CSV.read(input_file2, DataFrame)
    global x2 = DF2[:, Symbol(k)]
    global y2 = DF2[:, Symbol("f($k)")] ./ DF2[:, Symbol("Population Size") ]
    if maximum(x) <= x_thresh
        Plots.plot!(p1[k], x2, y2; seriestype=:scatter, 
                                        marker = marker[2], markercolor = mc[2], markerstrokecolor = mc[2],
                                        alpha=0.8,
                                        label="$label1 (Binomial)"
                                    )
    else
        Plots.plot!(p1[k], x2, y2; linestyle=:dash, color = mc[2], label="$label1 (Binomial)")
    end
    
    PyPlot.figure(k)
    a = (y.>0)
    a2 = (y2.>0)
    global p2[k] = Plots.plot( x[a], y[a];  label="$label1 (Uniform)", yscale=:log10)
    Plots.plot!(p2[k], x2[a2], y2[a2]; label="$label1 (Binomial)")
    Plots.plot!(p2[k]; xlabel="k", ylabel="log PMF")

    PyPlot.figure(k)
    cdf =  cumsum(y)
    cdf2 = cumsum(y2)
    global p3[k] = Plots.plot( x[2:end], 1.0 .- cdf[1:end-1];  label="$label1 (Uniform)", xscale=:log10, yscale=:log10)
    Plots.plot!(p3[k], x2[2:end], 1.0 .- cdf2[1:end-1];  label="$label1 (Binomial)")
    Plots.plot!(p3[k]; xlabel="log k", ylabel="log CCDF")
    Plots.plot!(p3[k]; legend=:bottomleft)

    if k=="Generation"
        epsilon = 0.1
        if maximum(x2) < maximum(x)
            xd = x
        else
            xd = x2
        end
        global geometric_appox = λ .* (1 .- α) .* α .^ xd
        pn = Int64.(round.(n .* geometric_appox  ))
        a = (geometric_appox .< epsilon)
        Plots.plot!(p1[k], xd[a], geometric_appox[a]; label="Geometric Approx", linestyle=:dash, )
    
        global pmf = generation_pmf.( xd; α=α, λ=λ)
        pm = Int64.(round.(n .* pmf))
        Plots.plot!(p1[k], xd, pmf; label="Predicted", linestyle=:dot, lw=2*lw)

        if maximum(xd) <= 15
            Plots.plot!(p1[k]; xticks = 0:2:maximum(xd))
        end
    end
    if k=="Log2_Denominator"
        f(x) = (2/x) * (1 - exp(-x) ) - exp(-x)
        # Plots.plot!(p1[k],  x2, exp(-λ).* λ.^x2 ./ factorial.(x2); label="Poisson parents")
        p_geom1 = 1 - y[1]
        p_geom2 = 1 - y2[1]
        Plots.plot!(p1[k],  x2, (1 - p_geom1).* p_geom1.^x2; linestyle=:dot, label="Geoemtric approx (Uniform)")
        Plots.plot!(p1[k],  x2, (1 - p_geom2).* p_geom2.^x2; linestyle=:dot, label="Geoemtric approx (Binomial)")
        Plots.plot!(p1[k],  [0.0, maximum(x2)], [f(λ), f(λ)]; linestyle=:dash, label=L"Lower bound for $k=0$")
    end

    global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_%s_uni_v_bin.pdf", distributions, λ,  α , n, g_0, PopulationSize, k))
    Plots.savefig(p1[k], plot_file)
    println("  output $plot_file")

    global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_%s_uni_v_bin_logy.pdf", distributions, λ,  α , n, g_0, PopulationSize, k))
    Plots.savefig(p2[k], plot_file)
    println("  output $plot_file")

    global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_%s_uni_v_bin_loglog.pdf", distributions, λ,  α , n, g_0, PopulationSize, k))
    Plots.savefig(p3[k], plot_file)
    println("  output $plot_file")

end 

p1["Generation"]

