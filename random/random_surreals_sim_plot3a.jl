# do some plots from the output of random_surreals_sim.jl
#  this one shows the evolutions through the iterations of average (and maybe max) stats
#    compare the impact of the size of the population
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

α = 0.8; λ = 3.5
# α = 0.4; λ = 1.5

g_0 = 1
Nseeds = 30

distributions = "P-U"
p2 = Vector{Any}(undef,3)

for (j,measure) in enumerate(["Generation", "Nodes", "Edges"])
    p2[j]  = Plots.plot( [0], [0]; xlabel="Iteration", ylabel="$measure",
                         legend=:outerright, label="")

    for n in Ns
        stat = "Mean"
        global p2, DF, filename, x, y_symbol_1, y_symbol_2, stop_state, m, m_final_states, k, ss, rel_difference, PopulationSize, base_file_name

        PopulationSize = n * Nseeds
        base_file_name = joinpath( input_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%05d_%02d_%07d_proc", 
                                                        distributions, λ,  α , n, g_0, PopulationSize))

        filename = "$(base_file_name)_$(measure)_$(stat).csv"
        println("  reading $measure $stat $filename")
        DF = CSV.read(filename, DataFrame)
            
        global m = size(DF,1)
        x = 0 : m-1
        y_symbol_1 = Symbol("Mean $stat $measure")
        y_symbol_2 = Symbol("Std $stat $measure")
        # p1[i,j] = Plots.plot( x, DF[:, y_symbol_1] ; label="", xlabel="iteration", ylabel="$stat -- $measure")

        Plots.plot!( p2[j], x, DF[:, y_symbol_1] ; label="n=$n")

        # stop_state = mean( DF[m-m_final_states+1 : m, y_symbol_1] )
        # ss = repeat( [stop_state], length(x) )
        # Plots.plot!( p2[j], x, ss ; label="", linestyle=:dot, linewidth=lw)

        # rel_difference = abs.( ss .- DF[:, y_symbol_1]) ./ stop_state
        # k = findfirst( rel_difference .< ϵ )
        # if !isnothing(k)
        #     Plots.plot!( p2[j], [k-1; k-1], [0; stop_state] ; label="", linestyle=:dot, linewidth=2*lw)
        # end

        
    end
    Plots.plot!( p2[j]; xlim=(0,m-1), ylim=(0,Inf) )

    global plot_file = joinpath( plot_dir, @sprintf("random_surreal_sim_%s_%0.2f_%0.2f_%02d_proc_3a_%s.pdf", 
                                                     distributions, λ, α, g_0, measure))
    Plots.savefig(p2[j], plot_file)

end

p2[1]

# do a version b that compares g_0
# do a version that plots Max and Mean on samples
